# =============================================================================
# vrmol_final_analysis.R
# Bachelor thesis: anxiety as a moderator of presence, memory and
#                  cybersickness in VRMoL (Method of Loci in VR)
#
# RQ1: presence (IPQ) -> memory (KoL free recall), moderated by anxiety
#      (DASS-21 anxiety), controlling for baseline memory (RAVLT) and
#      immersion group.
# RQ2: presence -> cybersickness (SSQ total change, post - pre), same
#      moderator and covariates.
#
# Parts:
#   1. Read the raw data, clean it and compute the scale scores
#   2. Descriptives, plots and regression assumption checks
#   3. Hierarchical moderated regressions, robustness checks, sensitivity reruns
#   4. Self-check against the step-by-step analysis + sessionInfo()
#
# All decisions (exclusions, missing data, recoding, scoring, sensitivity sets)
# are written in this script. Participants are identified by study ID only,
# with reason categories.
#
# Input : data/raw/kingdom_of_loci_raw.csv  
# Output: output/final/  
# Run from the project folder:
#   Rscript.exe scripts/vrmol_final_analysis.R
# =============================================================================


# =============================================================================
# 0. SETUP
# =============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(tibble)
  library(ggplot2)
})
# Called with :: (not loaded): psych, lmtest, car, boot, broom, interactions

raw_file <- "data/raw/kingdom_of_loci_raw.csv"
out_dir  <- "output/final"
out_txt  <- file.path(out_dir, "results.txt")

seed   <- 20261012   # fixed seed: the bootstrap gives the same numbers every run 
n_boot <- 5000       # bootstrap resamples
alpha  <- .05        # all tests two-tailed

if (!file.exists(raw_file)) stop("Raw file not found: ", raw_file, ". Run from the project folder.")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

options(width = 200)
sink(out_txt, split = TRUE)
cat("vrmol_final_analysis.R run on", format(Sys.time(), "%Y-%m-%d %H:%M"), "\n")

section <- function(title) {
  cat("\n", strrep("=", 78), "\n", title, "\n", strrep("=", 78), "\n", sep = "")
}
show <- function(x) print(as_tibble(x), n = Inf, width = Inf)
stop_here <- function(...) {
  sink()
  stop(..., call. = FALSE)
}
out_file <- function(name) file.path(out_dir, name)


# =============================================================================
# 0b. DECISIONS ABOUT PARTICIPANTS (IDs + reason categories only)
# =============================================================================

# Did not consent: removed from all analyses
no_consent_ids <- c("AU26EA", "NA08AN")

# Excluded from one research question only (NA = not excluded)
exclusions <- tribble(
  ~id,      ~rq,   ~reason,
  "HY02MJ", "RQ1", "free recall done without VR",
  "AR1506", "RQ2", "sick before VR (unrelated cause)"
)

# Score set to missing: session stopped during cued recall (secondary outcome)
cued_invalid_ids <- "OS28TS"

# Robustness check C: |standardized residual| > 3 in the RQ2 main model
rq2_outlier_ids <- c("AA01PE", "YI16TV")

# Sensitivity sets (reason category per ID), from the review of the session
# comments. Each set is left out once in a rerun of the final models (part 3).
#   protocol               protocol deviation during VR or learning
#   comprehension          questionnaire comprehension or language
#   characteristics        participant characteristics (prior MoL knowledge,
#                          prior RAVLT exposure, health or substance use)
#   ravlt_delay_only       RAVLT delay interval not standard (only used in the
#                          RAVLT-delay check)
#   la19rw_invalid_session session marked not valid by the experimenter (RQ2)
sensitivity_sets <- bind_rows(
  tibble(set = "protocol",
         id  = c("AE27LC", "AI08AE", "AO21DR", "AR20ZG", "AR30PV", "EA19MW", "EH20LZ", "EV30EX",
                 "FY24AA", "HY02MJ", "IE22EM", "IT27LH", "IW31SE", "LE12JW", "NN15LX", "OI19BS",
                 "OK19BD", "OS26AB", "RA07BG", "RK26HL", "RN26XV", "UA11LA", "US24TB")),
  tibble(set = "comprehension",
         id  = c("AA02BA", "AA02NM", "EH20LZ", "HH23RK", "HW04VM", "LO10HP")),
  tibble(set = "characteristics",
         id  = c("AA02BA", "AA13EM", "AR1506", "AR30PV", "EE19KG", "EI03RW", "IW31SE", "LA13HH", "LU06LX")),
  tibble(set = "ravlt_delay_only",       id = c("AL23JR", "IE22EM", "LA23AA")),
  tibble(set = "la19rw_invalid_session", id = "LA19RW")
)


# =============================================================================
# PART 1. READ, CLEAN AND SCORE
# =============================================================================

# ---- 1.1 Read the raw file and repair split records --------------------------
section("1.1 Read the raw data")

# Some answers in the raw export contain line breaks without quotes, which
# splits one participant over several lines. Lines are glued back together
# (with " / ") until a record has as many columns as the header.

# Number of columns in one line of CSV text (commas inside quotes not counted);
# NA if a quoted text is still open, i.e. the record continues on the next line
count_fields <- function(txt) {
  if (str_count(txt, fixed('"')) %% 2 == 1) return(NA_integer_)
  str_count(str_remove_all(txt, '"[^"]*"'), fixed(",")) + 1L
}

lines    <- read_lines(raw_file)
lines[1] <- str_remove(lines[1], "^﻿")   # remove invisible BOM character
header   <- lines[1]
n_cols   <- count_fields(header)

rec_text <- character()
buffer   <- NULL
for (i in 2:length(lines)) {
  buffer <- if (is.null(buffer)) lines[i] else paste(buffer, lines[i], sep = " / ")
  n <- count_fields(buffer)
  if (!is.na(n) && n >= n_cols) {
    rec_text <- c(rec_text, buffer)
    buffer   <- NULL
  }
}
if (!is.null(buffer)) stop_here("The last record is incomplete.")
n_fields <- map_int(rec_text, count_fields)
if (any(is.na(n_fields) | n_fields != n_cols)) stop_here("Records with a wrong number of columns.")

# Everything as text; missing codes are handled in 1.4
raw_all <- read_csv(I(c(header, rec_text)), col_types = cols(.default = col_character()),
                    na = character(), trim_ws = TRUE, progress = FALSE)
if (nrow(problems(raw_all)) > 0) stop_here("read_csv reported parsing problems.")

# Participant rows only (drops the footer row)
raw <- raw_all %>% filter(str_detect(IDSPSS, "^[A-Z]{2}[0-9A-Z]{2}[A-Z0-9]{2}$"))
if (anyDuplicated(raw$IDSPSS) > 0) stop_here("Duplicated participant IDs.")
cat(length(lines), "lines ->", nrow(raw_all), "records x", n_cols, "columns ->",
    nrow(raw), "participants (expected 109)\n")
if (nrow(raw) != 109) stop_here("Expected 109 participants.")


# ---- 1.2 Remove participants who did not consent -----------------------------
section("1.2 Consent")

no_consent <- raw %>% filter(str_detect(InfoConsent, "DO NOT AGREE")) %>% pull(IDSPSS)
if (!setequal(no_consent, no_consent_ids)) stop_here("Non-consenters differ from the decision list.")
consented <- raw %>% filter(!IDSPSS %in% no_consent_ids)
cat("Did not consent (removed):", paste(no_consent, collapse = ", "),
    "| participants left:", nrow(consented), "(expected 107)\n")


# ---- 1.3 Keep only the needed variables, with clear names --------------------
section("1.3 Select variables")

anx_items <- c(2, 4, 7, 9, 15, 19, 20)   # DASS-21 anxiety items

# new name = column in the raw file. Questionnaires use the file's own
# numeric _Recoded columns (no recoding twice).
keep <- c(
  id          = "IDSPSS",
  age         = "Age",
  sex         = "Sex",                   # 0 = Female, 1 = Male (.sav value labels)
  education   = "Educationlevel",        # highest completed education (answer options)
  in_education = "Currenteducation",     # currently in education: Yes / No
  nationality = "Nationality",           # free text; grouped in 2.2
  language    = "Primarylanguage",       # 0 = Dutch, 1 = English, 3 = Other; several allowed
  group_s1    = "S1Condition",
  group       = "S2Condition",           # immersion group (covariate)
  set_names(paste0("RAVLT", 1:5), paste0("ravlt_", 1:5)),
  ravlt_delay = "RAVLTDelay",
  set_names(paste0("DASS_", anx_items, "_Recoded"), paste0("dass_", anx_items)),
  set_names(paste0("IPQ", 1:14, "_Recoded"), paste0("ipq_", 1:14)),
  set_names(paste0("SSQPre_", 1:16, "_Recoded"), paste0("ssq_pre_", 1:16)),
  set_names(paste0("SSQPOST_", 1:16, "_Recoded"), paste0("ssq_post_", 1:16)),
  fms_1 = "FMSickness1.1", fms_2 = "FMSickness1.2", fms_3 = "FMSickness1.3",
  fms_4 = "FMSickness2.1", fms_5 = "FMSickness2.2", fms_6 = "FMSickness3.2",
  fms_7 = "FMSickness3.3",
  kol_free_r1   = "KOLfreeR1",           # RQ1 outcome
  kol_cued_r    = "KOLcuedR",            # secondary outcome
  kol_cued_mp   = "KOLcuedMP",
  kol_cued_pass = "KOLcuedPASS",
  kol_cued_intr = "KOLcuedINTR"
)
missing_cols <- setdiff(keep, names(consented))
if (length(missing_cols) > 0) stop_here("Columns not found: ", paste(missing_cols, collapse = ", "))

selected <- consented %>%
  select(all_of(keep)) %>%
  mutate(across(everything(), str_trim))
cat("Kept", ncol(selected), "of", ncol(consented), "columns.\n")


# ---- 1.4 Missing codes -> NA -------------------------------------------------
section("1.4 Missing codes (999, 99, empty) -> NA")

missing_codes <- c("999", "99", "")
show(selected %>%
       select(-id) %>%
       pivot_longer(everything(), names_to = "variable", values_to = "value") %>%
       group_by(variable) %>%
       summarise(n_empty = sum(value == ""), n_999 = sum(value == "999"),
                 n_99 = sum(value == "99"), .groups = "drop") %>%
       filter(n_empty + n_999 + n_99 > 0))

cleaned <- selected %>%
  mutate(across(-id, ~ if_else(.x %in% missing_codes, NA_character_, .x)))


# ---- 1.5 Numbers, rescaling and the group factor -----------------------------
section("1.5 Convert to numbers; group factor")

to_number <- function(x, var) {
  num <- suppressWarnings(as.numeric(x))
  if (any(!is.na(x) & is.na(num))) stop_here("Non-numeric value(s) in ", var)
  num
}

numeric_vars <- c("age", paste0("ravlt_", 1:5), "ravlt_delay", paste0("fms_", 1:7),
                  str_subset(names(cleaned), "^(kol|dass|ssq|ipq)_"))

converted <- cleaned %>%
  mutate(
    across(all_of(numeric_vars), ~ to_number(.x, cur_column())),
    # DASS and SSQ: file codes 1-4 minus 1 -> standard 0-3 (order unchanged)
    across(starts_with(c("dass_", "ssq_")), ~ .x - 1),
    # IPQ: file codes -3..+3 used as they are (items 2, 11, 13 reversed in the file)
    # Immersion group: factor, low (code 3) = reference category
    group    = factor(group,    levels = c("3", "1", "2"), labels = c("low", "high", "moderate")),
    group_s1 = factor(group_s1, levels = c("3", "1", "2"), labels = c("low", "high", "moderate")),
    # Sex: factor with Female (code 0) as reference category
    sex      = factor(sex, levels = c("0", "1"), labels = c("Female", "Male"))
  )
show(count(converted, group))
cat("Session 1 group differs from session 2 group for:",
    paste(converted$id[which(converted$group_s1 != converted$group)], collapse = ", "), "\n")


# ---- 1.6 Participant decisions -----------------------------------------------
section("1.6 Participant decisions")

if (!all(c(exclusions$id, cued_invalid_ids) %in% converted$id)) stop_here("Decision IDs not found.")

decided <- converted %>%
  mutate(kol_cued_r = if_else(id %in% cued_invalid_ids, NA_real_, kol_cued_r)) %>%
  left_join(exclusions %>% filter(rq == "RQ1") %>% select(id, excl_rq1 = reason), by = "id") %>%
  left_join(exclusions %>% filter(rq == "RQ2") %>% select(id, excl_rq2 = reason), by = "id")
show(exclusions)
cat("Cued recall set to missing for:", cued_invalid_ids, "\n")


# ---- 1.7 Checks: unique IDs and valid ranges ---------------------------------
section("1.7 Checks: unique IDs and valid item ranges")

if (anyDuplicated(decided$id) > 0) stop_here("Duplicated IDs.")

ranges <- tribble(
  ~pattern,                     ~lo, ~hi,
  "^age$",                       18,  80,
  "^ravlt_([0-9]|delay)$",        0,  15,
  "^dass_[0-9]+$",                0,   3,
  "^ssq_(pre|post)_[0-9]+$",      0,   3,
  "^ipq_[0-9]+$",                -3,   3,
  "^fms_[0-9]$",                  0,  20,
  "^kol_",                        0,  27
)
out_of_range <- decided %>%
  select(id, matches(paste(ranges$pattern, collapse = "|"))) %>%
  pivot_longer(-id, names_to = "variable", values_to = "value") %>%
  filter(!is.na(value)) %>%
  mutate(row = map_int(variable, ~ which(str_detect(.x, ranges$pattern))[1]),
         lo = ranges$lo[row], hi = ranges$hi[row]) %>%
  filter(value < lo | value > hi | value != round(value))
cat("Item values outside their range or not whole:", nrow(out_of_range), "\n")
if (nrow(out_of_range) > 0) { show(out_of_range); stop_here("Values outside the valid range.") }


# ---- 1.8 Scale scores --------------------------------------------------------
section("1.8 Scale scores")

# A score is only computed when ALL its items are present (rowSums / rowMeans
# give NA otherwise); the FMS peak uses the available ratings.
ravlt_items <- paste0("ravlt_", 1:5)
dass_items  <- paste0("dass_", anx_items)
ipq_items   <- paste0("ipq_", 1:14)
fms_items   <- paste0("fms_", 1:7)

# SSQ subscales (Kennedy et al., 1993), items in the standard order
ssq_n <- c(1, 6, 7, 8, 9, 15, 16)      # Nausea
ssq_o <- c(1, 2, 3, 4, 5, 9, 11)       # Oculomotor
ssq_d <- c(5, 8, 10, 11, 12, 13, 14)   # Disorientation

ssq_scores <- function(df, time) {
  items <- function(nums) df[paste0("ssq_", time, "_", nums)]
  n_raw <- rowSums(items(ssq_n))
  o_raw <- rowSums(items(ssq_o))
  d_raw <- rowSums(items(ssq_d))
  tibble(n = n_raw * 9.54, o = o_raw * 7.58, d = d_raw * 13.92,
         total = (n_raw + o_raw + d_raw) * 3.74) %>%
    rename_with(~ paste0("ssq_", time, "_", .x))
}

d <- decided %>%
  mutate(
    ravlt_total = rowSums(across(all_of(ravlt_items))),        # 0-75
    dass_anx    = rowSums(across(all_of(dass_items))) * 2,     # 0-42 (sum x 2)
    ipq_mean    = rowMeans(across(all_of(ipq_items))),         # -3 to +3
    fms_n       = rowSums(!is.na(across(all_of(fms_items)))),
    fms_peak    = if_else(fms_n > 0, pmax(!!!syms(fms_items), na.rm = TRUE), NA_real_)
  ) %>%
  bind_cols(ssq_scores(decided, "pre"), ssq_scores(decided, "post")) %>%
  mutate(ssq_change = ssq_post_total - ssq_pre_total)

ssq_max <- 3 * 21 * 3.74
limits <- tribble(
  ~score,           ~lo,       ~hi,
  "ravlt_total",      0,        75,
  "dass_anx",         0,        42,
  "ipq_mean",        -3,         3,
  "ssq_pre_total",    0,   ssq_max,
  "ssq_post_total",   0,   ssq_max,
  "ssq_change", -ssq_max,  ssq_max,
  "fms_peak",         0,        20
)
bad_scores <- d %>%
  select(id, all_of(limits$score)) %>%
  pivot_longer(-id, names_to = "score", values_to = "value") %>%
  left_join(limits, by = "score") %>%
  filter(!is.na(value), value < lo | value > hi)
cat("Scores outside their possible range:", nrow(bad_scores), "\n")
if (nrow(bad_scores) > 0) stop_here("Scores outside the possible range.")

score_vars <- c("ravlt_total", "dass_anx", "ipq_mean", "kol_free_r1", "kol_cued_r",
                "ssq_pre_total", "ssq_post_total", "ssq_change", "fms_peak")
show(d %>%
       select(all_of(score_vars)) %>%
       pivot_longer(everything(), names_to = "score", values_to = "value") %>%
       group_by(score = factor(score, levels = score_vars)) %>%
       summarise(n = sum(!is.na(value)), M = mean(value, na.rm = TRUE), SD = sd(value, na.rm = TRUE),
                 min = min(value, na.rm = TRUE), max = max(value, na.rm = TRUE), .groups = "drop") %>%
       mutate(across(c(M, SD, min, max), ~ round(.x, 2))))


# ---- 1.9 Reliability: Cronbach's alpha ---------------------------------------
section("1.9 Reliability (Cronbach's alpha, complete cases per scale)")

run_alpha <- function(scale, items) {
  item_data <- d %>% select(all_of(items)) %>% drop_na()
  a <- suppressWarnings(psych::alpha(as.data.frame(item_data), check.keys = FALSE))
  tibble(scale = scale, n = nrow(item_data), items = length(items),
         alpha = round(a$total$raw_alpha, 3))
}
alpha_table <- bind_rows(
  run_alpha("IPQ (14 items)",         ipq_items),
  run_alpha("DASS anxiety (7 items)", dass_items),
  run_alpha("SSQ pre (16 items)",     paste0("ssq_pre_", 1:16)),
  run_alpha("SSQ post (16 items)",    paste0("ssq_post_", 1:16))
)
show(alpha_table)
write_excel_csv(alpha_table, out_file("1_reliability.csv"))
cat("Scored dataset kept in memory only (", nrow(d), "x", ncol(d), "); no participant-level file is saved.\n")


# =============================================================================
# PART 2. DESCRIPTIVES, PLOTS AND ASSUMPTION CHECKS
# =============================================================================

# ---- 2.1 Analysis samples ----------------------------------------------------
section("2.1 Analysis samples (complete cases per model + decision exclusions)")

needed_both <- c("ipq_mean", "dass_anx", "ravlt_total", "group")
d <- d %>%
  mutate(in_rq1 = if_all(all_of(c(needed_both, "kol_free_r1")), ~ !is.na(.x)) & is.na(excl_rq1),
         in_rq2 = if_all(all_of(c(needed_both, "ssq_change")),  ~ !is.na(.x)) & is.na(excl_rq2))
rq1 <- d %>% filter(in_rq1)
rq2 <- d %>% filter(in_rq2)
cat("RQ1 N =", nrow(rq1), "| RQ2 N =", nrow(rq2), "| in both:", sum(d$in_rq1 & d$in_rq2), "\n")

# Sample flow: participants missing each variable (can overlap) and excluded
flow <- function(outcome, excl_col, in_col, label) {
  needed <- c(needed_both, outcome)
  tibble(sample = label,
         reason = c(paste("missing", needed), "excluded by decision", "lost in total", "final N"),
         n      = c(map_int(needed, ~ sum(is.na(d[[.x]]))), sum(!is.na(d[[excl_col]])),
                    sum(!d[[in_col]]), sum(d[[in_col]])))
}
show(bind_rows(flow("kol_free_r1", "excl_rq1", "in_rq1", "RQ1"),
               flow("ssq_change",  "excl_rq2", "in_rq2", "RQ2")))


# ---- 2.2 Demographics --------------------------------------------------------
section("2.2 Demographics")

samples <- list("All consenters" = d, "RQ1" = rq1, "RQ2" = rq2)
demo <- function(s) {
  summarise(s, N = n(), age_n = sum(!is.na(age)),
            age_M = round(mean(age, na.rm = TRUE), 2), age_SD = round(sd(age, na.rm = TRUE), 2),
            age_min = min(age, na.rm = TRUE), age_max = max(age, na.rm = TRUE),
            sex_female = sum(sex == "Female", na.rm = TRUE), sex_male = sum(sex == "Male", na.rm = TRUE),
            sex_missing = sum(is.na(sex)), .groups = "drop")
}
show(imap_dfr(samples, ~ demo(.x) %>% mutate(sample = .y)) %>% relocate(sample))

cat("\nAge and sex per immersion group (descriptive only, no tests):\n")
group_demographics <- imap_dfr(list(RQ1 = rq1, RQ2 = rq2),
                               ~ demo(group_by(.x, group)) %>% mutate(sample = .y)) %>%
  relocate(sample)
show(group_demographics)
write_excel_csv(group_demographics, out_file("2_group_demographics.csv"))

# Background for the Participants section: education, nationality, language.
# Nationality is free text: spelling variants of the same nationality are
# merged (lower case; country or Dutch-language names -> nationality), then
# only Dutch / other / missing and the number of distinct other nationalities
# are reported.
nationality_merge <- c(nederlands = "dutch", nederlandse = "dutch", spain = "spanish",
                       taiwan = "taiwanese", belgisch = "belgian")
background <- function(s) {
  s <- s %>%
    mutate(nat   = str_to_lower(nationality),
           nat   = coalesce(unname(nationality_merge[nat]), nat),
           # Primary language: code 0 or the text "Dutch" (multiple answers)
           dutch = str_detect(language, "^0$|Dutch"),
           only_dutch = language %in% c("0", "Dutch"))
  other_nat <- s$nat[!is.na(s$nat) & s$nat != "dutch"]
  bind_rows(
    s %>% count(category = coalesce(education, "missing")) %>% mutate(measure = "education level"),
    s %>% count(category = coalesce(in_education, "missing")) %>% mutate(measure = "currently in education"),
    tibble(measure  = "nationality",
           category = c("Dutch", "other", "missing", "distinct other nationalities"),
           n = c(sum(s$nat == "dutch", na.rm = TRUE), length(other_nat), sum(is.na(s$nat)),
                 n_distinct(other_nat))),
    tibble(measure  = "primary language",
           category = c("Dutch (alone or with another)", "Dutch only", "no Dutch", "missing"),
           n = c(sum(s$dutch, na.rm = TRUE), sum(s$only_dutch), sum(!s$dutch, na.rm = TRUE),
                 sum(is.na(s$language))))
  ) %>%
    select(measure, category, n)
}
background_table <- imap_dfr(samples, ~ background(.x) %>% mutate(sample = .y)) %>%
  pivot_wider(names_from = sample, values_from = n, values_fill = 0L)
cat("\nEducation, nationality and primary language (n per sample):\n")
show(background_table)
write_excel_csv(background_table, out_file("2_background.csv"))


# ---- 2.3 Descriptives of the model variables ---------------------------------
section("2.3 Descriptives of the model variables (overall and per group)")

desc_vars <- c(
  ipq_mean       = "Presence (IPQ mean)",
  dass_anx       = "Anxiety (DASS-21 x 2)",
  ravlt_total    = "RAVLT total learning",
  kol_free_r1    = "KoL free recall (RQ1 outcome)",
  ssq_change     = "SSQ change (RQ2 outcome)",
  kol_cued_r     = "KoL cued recall (secondary)",
  ssq_pre_total  = "SSQ pre total",
  ssq_post_total = "SSQ post total",
  fms_peak       = "FMS peak (secondary)"
)
describe_one <- function(x) {
  x <- x[!is.na(x)]
  tibble(n = length(x), M = mean(x), SD = sd(x), min = min(x), max = max(x),
         skew = psych::skew(x), kurtosis = psych::kurtosi(x))
}
describe_sample <- function(s, name) {
  bind_rows(s %>% mutate(group = "all"), s %>% mutate(group = as.character(group))) %>%
    select(group, all_of(names(desc_vars))) %>%
    pivot_longer(-group, names_to = "variable", values_to = "value") %>%
    group_by(group, variable) %>%
    summarise(describe_one(value), .groups = "drop") %>%
    mutate(sample   = name,
           variable = factor(variable, levels = names(desc_vars), labels = desc_vars),
           group    = factor(group, levels = c("all", "low", "high", "moderate"))) %>%
    arrange(variable, group) %>%
    relocate(sample, variable, group)
}
desc_table <- imap_dfr(samples, describe_sample) %>%
  mutate(across(c(M, SD, min, max, skew, kurtosis), ~ round(.x, 2)))
for (name in names(samples)) {
  cat("\n---", name, "---\n")
  show(desc_table %>% filter(sample == name) %>% select(-sample))
}
write_excel_csv(desc_table, out_file("2_descriptives.csv"))


# ---- 2.4 Correlations --------------------------------------------------------
section("2.4 Pearson correlations")

cor_out <- function(s, vars, file, label) {
  r <- as_tibble(round(cor(s[vars]), 2), rownames = "variable")
  cat("\n", label, " (n = ", sum(complete.cases(s[vars])), "):\n", sep = "")
  show(r)
  write_excel_csv(r, out_file(file))
}
cor_out(rq1, c("ipq_mean", "dass_anx", "ravlt_total", "kol_free_r1"), "2_correlations_rq1.csv", "RQ1 sample")
cor_out(rq2, c("ipq_mean", "dass_anx", "ravlt_total", "ssq_change"),  "2_correlations_rq2.csv", "RQ2 sample")
both <- d %>% filter(in_rq1, in_rq2)
cor_out(both, c("ipq_mean", "dass_anx", "ravlt_total", "kol_free_r1", "ssq_change"),
        "2_correlations_all5.csv", "All five variables, participants in both samples")

# Exploratory: do the two cybersickness measures agree (FMS peak during VR vs
# SSQ change)? Spearman added because SSQ change is right-skewed.
fms_ssq <- rq2 %>% filter(!is.na(fms_peak))
fms_pearson  <- cor.test(fms_ssq$fms_peak, fms_ssq$ssq_change, method = "pearson")
fms_spearman <- cor.test(fms_ssq$fms_peak, fms_ssq$ssq_change, method = "spearman", exact = FALSE)
fms_table <- tibble(
  sample = "RQ2 (with FMS ratings)", n = nrow(fms_ssq),
  pearson_r = fms_pearson$estimate, ci_low = fms_pearson$conf.int[1], ci_high = fms_pearson$conf.int[2],
  p_pearson = fms_pearson$p.value, spearman_rho = fms_spearman$estimate, p_spearman = fms_spearman$p.value
)
cat("\nFMS peak x SSQ change (exploratory):\n")
show(fms_table %>% mutate(across(c(pearson_r, ci_low, ci_high, spearman_rho), ~ round(.x, 3)),
                          across(starts_with("p_"), ~ signif(.x, 3))))
write_excel_csv(fms_table, out_file("2_fms_ssq_correlation.csv"))


# ---- 2.5 Descriptive plots ---------------------------------------------------
section("2.5 Descriptive plots")

# Histograms of every model variable (all consenters with a value)
p_hist <- d %>%
  select(all_of(names(desc_vars))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "value") %>%
  filter(!is.na(value)) %>%
  mutate(variable = factor(variable, levels = names(desc_vars), labels = desc_vars)) %>%
  ggplot(aes(x = value)) +
  geom_histogram(bins = 15, fill = "grey60", colour = "white") +
  facet_wrap(~ variable, scales = "free", ncol = 3) +
  labs(title = "Distributions of the model variables (consenting participants)",
       x = NULL, y = "Number of participants") +
  theme_minimal()
ggsave(out_file("2_histograms.png"), p_hist, width = 10, height = 8, dpi = 300)

# Boxplots per immersion group (main model variables)
p_box <- d %>%
  select(group, ipq_mean, dass_anx, kol_free_r1, ssq_change) %>%
  pivot_longer(-group, names_to = "variable", values_to = "value") %>%
  filter(!is.na(value)) %>%
  mutate(variable = droplevels(factor(variable, levels = names(desc_vars), labels = desc_vars))) %>%
  ggplot(aes(x = group, y = value)) +
  geom_boxplot(outlier.shape = NA, fill = "grey90") +
  geom_jitter(width = 0.15, alpha = 0.5, size = 1) +
  facet_wrap(~ variable, scales = "free_y", ncol = 2) +
  labs(title = "Main model variables per immersion group", x = "Immersion group", y = NULL) +
  theme_minimal()
ggsave(out_file("2_boxplots_group.png"), p_box, width = 9, height = 7, dpi = 300)

# Presence vs outcome per analysis sample (raw relationship, not the model)
p_scatter <- bind_rows(
  rq1 %>% transmute(panel = "RQ1: presence vs KoL free recall", ipq_mean, outcome = kol_free_r1),
  rq2 %>% transmute(panel = "RQ2: presence vs SSQ change", ipq_mean, outcome = ssq_change)
) %>%
  ggplot(aes(x = ipq_mean, y = outcome)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = "black") +
  facet_wrap(~ panel, scales = "free_y") +
  labs(title = "Presence and the two outcomes (raw relationship)",
       x = "Presence (IPQ mean)", y = "Outcome") +
  theme_minimal()
ggsave(out_file("2_scatter_presence.png"), p_scatter, width = 10, height = 4.5, dpi = 300)
cat("Saved: 2_histograms.png, 2_boxplots_group.png, 2_scatter_presence.png\n")


# ---- 2.6 Centring and the full models ----------------------------------------
section("2.6 Centring (presence and anxiety, within each RQ sample) and full models")

# Subtract the sample mean, so 0 = the average participant of that analysis.
# RAVLT stays a raw score; group stays a factor (reference = low).
centre <- function(s) {
  s %>% mutate(presence_c = ipq_mean - mean(ipq_mean),
               anxiety_c  = dass_anx - mean(dass_anx))
}
rq1 <- centre(rq1) %>% mutate(y = kol_free_r1)
rq2 <- centre(rq2) %>% mutate(y = ssq_change)
show(tibble(sample = c("RQ1", "RQ2"),
            presence_mean = c(mean(rq1$ipq_mean), mean(rq2$ipq_mean)),
            anxiety_mean  = c(mean(rq1$dass_anx), mean(rq2$dass_anx))))

f_step1 <- y ~ ravlt_total + group
f_step2 <- y ~ ravlt_total + group + presence_c + anxiety_c
f_step3 <- y ~ ravlt_total + group + presence_c + anxiety_c + presence_c:anxiety_c

full_models <- list(RQ1 = lm(f_step3, data = rq1), RQ2 = lm(f_step3, data = rq2))
rq_samples  <- list(RQ1 = rq1, RQ2 = rq2)


# ---- 2.7 Assumption checks ---------------------------------------------------
section("2.7 Assumption checks on the full models")

cat("Normality of residuals (Shapiro-Wilk):\n")
show(imap_dfr(full_models, function(m, name) {
  sw <- shapiro.test(residuals(m))
  tibble(model = name, W = round(sw$statistic, 3), p = signif(sw$p.value, 3),
         resid_skew = round(psych::skew(residuals(m)), 2),
         resid_kurtosis = round(psych::kurtosi(residuals(m)), 2))
}))

cat("\nHomoscedasticity (studentized Breusch-Pagan, Koenker 1981):\n")
show(imap_dfr(full_models, function(m, name) {
  bp <- lmtest::bptest(m, studentize = TRUE)
  tibble(model = name, BP = round(bp$statistic, 2), df = bp$parameter, p = signif(bp$p.value, 3))
}))

cat("\nMulticollinearity (VIF; GVIF for the group factor):\n")
for (name in names(full_models)) {
  v <- suppressMessages(suppressWarnings(car::vif(full_models[[name]])))
  v_tbl <- if (is.matrix(v)) as_tibble(v, rownames = "term") else tibble(term = names(v), VIF = v)
  cat(name, ":\n")
  show(v_tbl %>% mutate(across(where(is.numeric), ~ round(.x, 2))))
}

cat("\nOutliers and influence (|std. residual| > 3, Cook's D > 4/n, leverage > 2p/n; listed only):\n")
influence_table <- function(m, s, name) {
  n <- nobs(m); p <- length(coef(m))
  tibble(model = name, id = s$id, std_resid = rstandard(m),
         cooks_d = cooks.distance(m), leverage = hatvalues(m)) %>%
    mutate(flag_resid = abs(std_resid) > 3, flag_cook = cooks_d > 4 / n,
           flag_leverage = leverage > 2 * p / n)
}
infl <- imap_dfr(full_models, ~ influence_table(.x, rq_samples[[.y]], .y))
show(infl %>% group_by(model) %>%
       summarise(resid = sum(flag_resid), cook = sum(flag_cook), leverage = sum(flag_leverage),
                 max_cooks_d = round(max(cooks_d), 3)))
show(infl %>% filter(flag_resid | flag_cook | flag_leverage) %>%
       arrange(model, desc(cooks_d)) %>%
       mutate(across(c(std_resid, cooks_d, leverage), ~ round(.x, 3))))

# The RQ2 outliers used in robustness check C must be the |z| > 3 cases
found_outliers <- infl %>% filter(model == "RQ2", flag_resid) %>% pull(id)
if (!setequal(found_outliers, rq2_outlier_ids)) stop_here("RQ2 |z| > 3 cases differ from AA01PE, YI16TV.")
if (any(infl$model == "RQ1" & infl$flag_resid)) stop_here("RQ1 now has |z| > 3 cases.")


# ---- 2.8 Diagnostic plots ----------------------------------------------------
section("2.8 Diagnostic plots")

# Residuals vs fitted, Q-Q, scale-location, residuals vs leverage
for (name in names(full_models)) {
  png(out_file(paste0("2_diagnostics_", tolower(name), ".png")), width = 2400, height = 2400, res = 300)
  par(mfrow = c(2, 2))
  plot(full_models[[name]], which = c(1, 2, 3, 5), labels.id = rq_samples[[name]]$id)
  mtext(paste(name, "full model"), outer = TRUE, line = -1.5, font = 2)
  dev.off()
}

# Residuals vs each continuous predictor (linearity)
p_resid <- imap_dfr(rq_samples, function(s, name) {
  s %>%
    mutate(residual = residuals(full_models[[name]]), model = name) %>%
    select(model, residual, presence_c, anxiety_c, ravlt_total) %>%
    pivot_longer(c(presence_c, anxiety_c, ravlt_total), names_to = "predictor", values_to = "value")
}) %>%
  ggplot(aes(x = value, y = residual)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", formula = y ~ x, se = FALSE, colour = "black") +
  facet_grid(model ~ predictor, scales = "free") +
  labs(title = "Residuals vs each continuous predictor (full models)",
       x = "Predictor value", y = "Residual") +
  theme_minimal()
ggsave(out_file("2_resid_vs_predictors.png"), p_resid, width = 10, height = 6.5, dpi = 300)
cat("Saved: 2_diagnostics_rq1.png, 2_diagnostics_rq2.png, 2_resid_vs_predictors.png\n")


# =============================================================================
# PART 3. HYPOTHESIS TESTS, ROBUSTNESS AND SENSITIVITY
# =============================================================================

# Bias-corrected and accelerated (BCa) 95% interval for estimate number i of a
# bootstrap object. The bias correction (z0) adjusts for a shift of the
# bootstrap distribution away from the original estimate; the acceleration (a)
# adjusts for skewness (the SE changing with the size of the estimate).
bca_ci <- function(boot_out, i, level = .95) {
  ci <- suppressWarnings(boot::boot.ci(boot_out, conf = level, type = "bca", index = i))
  unname(ci$bca[1, 4:5])
}
outcome_label <- c(RQ1 = "KoL free recall", RQ2 = "SSQ change (post - pre)")


# ---- 3.1 Hierarchical regression ---------------------------------------------
section("3.1 Hierarchical regression: R2, delta R2, F change and f2 per step")

# Step 1: covariates; step 2: + presence + anxiety; step 3: + interaction
models <- map(rq_samples, function(s) {
  list(step1 = lm(f_step1, data = s), step2 = lm(f_step2, data = s), step3 = lm(f_step3, data = s))
})

step_table <- imap_dfr(models, function(m, name) {
  r2   <- map_dbl(m, ~ summary(.x)$r.squared)
  f1   <- summary(m$step1)$fstatistic             # step 1 vs the empty model
  comp <- anova(m$step1, m$step2, m$step3)         # steps 2 and 3 vs previous step
  tibble(
    sample   = name,
    step     = c("1: RAVLT + group", "2: + presence + anxiety", "3: + presence x anxiety"),
    R2       = r2,
    adj_R2   = map_dbl(m, ~ summary(.x)$adj.r.squared),
    delta_R2 = c(r2[1], diff(r2)),
    F_change = c(f1[["value"]], comp$F[2], comp$F[3]),
    df1      = c(f1[["numdf"]], comp$Df[2], comp$Df[3]),
    df2      = c(f1[["dendf"]], comp$Res.Df[2], comp$Res.Df[3]),
    p_change = c(pf(f1[["value"]], f1[["numdf"]], f1[["dendf"]], lower.tail = FALSE),
                 comp$`Pr(>F)`[2], comp$`Pr(>F)`[3]),
    f2       = delta_R2 / (1 - R2)                 # Cohen's f2 of the added block
  )
})
show(step_table %>% mutate(across(c(R2, adj_R2, delta_R2, f2), ~ round(.x, 3)),
                           F_change = round(F_change, 2), p_change = signif(p_change, 3)))
cat("f2 (Cohen, 1988): .02 = small, .15 = medium, .35 = large.\n")
write_excel_csv(step_table, out_file("3_steps.csv"))

# Sensitivity power analysis: smallest interaction f2 detectable with 80% power
# at alpha .05 with the fixed N (same model as G*Power "Linear multiple
# regression: Fixed model, R2 increase"; noncentral F with lambda = f2 * N).
power_f2 <- function(f2, n, n_tested, n_predictors, a = alpha) {
  df2 <- n - n_predictors - 1
  pf(qf(1 - a, n_tested, df2), n_tested, df2, ncp = f2 * n, lower.tail = FALSE)
}
power_table <- imap_dfr(models, function(m, name) {
  n <- nobs(m$step3)
  k <- length(coef(m$step3)) - 1                   # predictors in the final model
  min_f2 <- uniroot(function(f) power_f2(f, n, 1, k) - .80, c(1e-4, 2), tol = 1e-10)$root
  tibble(sample = name, N = n, tested = 1, predictors = k, alpha = alpha, power = .80,
         min_f2 = min_f2, observed_f2 = step_table$f2[step_table$sample == name][3])
})
cat("\nSensitivity power analysis (interaction, 80% power):\n")
show(power_table %>% mutate(across(c(min_f2, observed_f2), ~ round(.x, 3))))
write_excel_csv(power_table, out_file("3_power_sensitivity.csv"))


# ---- 3.2 Bootstrap of the final model ----------------------------------------
section("3.2 Bootstrap (5000 case resamples, seed 20261012)")

# One resample: refit step 3; return the coefficients and the simple slopes of
# presence at -1 SD, mean and +1 SD anxiety. Centring and anx_sd are fixed
# from the original sample (not recomputed per resample).
boot_stat <- function(s, idx, anx_sd) {
  b <- coef(lm(f_step3, data = s[idx, ]))
  slopes <- b[["presence_c"]] + b[["presence_c:anxiety_c"]] * c(-anx_sd, 0, anx_sd)
  c(b, slope_low = slopes[1], slope_mean = slopes[2], slope_high = slopes[3])
}
boots <- imap(rq_samples, function(s, name) {
  set.seed(seed)
  b <- boot::boot(s, boot_stat, R = n_boot, anx_sd = sd(s$anxiety_c))
  cat(name, ": resamples with a missing estimate:", sum(!complete.cases(b$t)), "\n")
  b
})


# ---- 3.3 Coefficients of the final model -------------------------------------
section("3.3 Final model (step 3): OLS and bootstrap BCa 95% CIs")

coef_table <- imap_dfr(models, function(m, name) {
  b <- boots[[name]]
  k <- length(coef(m$step3))
  bca <- map(seq_len(k), ~ bca_ci(b, .x))
  broom::tidy(m$step3, conf.int = TRUE) %>%
    mutate(sample   = name,
           bca_low  = map_dbl(bca, 1),
           bca_high = map_dbl(bca, 2)) %>%
    select(sample, term, b = estimate, SE = std.error, t = statistic, p = p.value,
           ols_low = conf.low, ols_high = conf.high, bca_low, bca_high)
})
for (name in names(models)) {
  cat("\n", name, " - ", outcome_label[[name]], " (n = ", nobs(models[[name]]$step3), ")\n", sep = "")
  show(coef_table %>% filter(sample == name) %>% select(-sample) %>%
         mutate(across(c(b, SE, t, ols_low, ols_high, bca_low, bca_high), ~ round(.x, 3)), p = signif(p, 3)))
}
cat("\nDecision rule: OLS p (two-tailed, alpha = .05); bootstrap BCa CIs as support.\n")
write_excel_csv(coef_table, out_file("3_coefficients.csv"))


# ---- 3.4 Simple slopes -------------------------------------------------------
section("3.4 Simple slopes of presence at -1 SD / mean / +1 SD anxiety")

# Interpreted only if the interaction (step 3) is significant. Johnson-Neyman
# is not reported (the main interactions are not significant; Nina, ACTION 67).
slope_table <- imap_dfr(models, function(m, name) {
  m3 <- m$step3
  p_interaction <- summary(m3)$coefficients["presence_c:anxiety_c", "Pr(>|t|)"]
  cat("\n---", name, "-", outcome_label[[name]], "--- interaction p =", signif(p_interaction, 3),
      if (p_interaction < alpha) "-> slopes interpreted\n" else "-> slopes shown, NOT interpreted\n")
  ss <- interactions::sim_slopes(m3, pred = presence_c, modx = anxiety_c,
                                 modx.values = "mean-plus-minus", centered = "none",
                                 johnson_neyman = FALSE)
  print(ss)
  k <- length(coef(m3))
  boot_ci <- map(1:3, ~ bca_ci(boots[[name]], k + .x))
  as_tibble(ss$slopes) %>%
    transmute(sample = name, anxiety_level = c("-1 SD", "mean", "+1 SD"),
              anxiety_c_value = .data[[names(ss$slopes)[1]]],
              slope = `Est.`, SE = `S.E.`, p = p, ols_low = `2.5%`, ols_high = `97.5%`,
              bca_low = map_dbl(boot_ci, 1), bca_high = map_dbl(boot_ci, 2),
              interaction_p = p_interaction, interpreted = p_interaction < alpha)
})
cat("\nSimple slopes (OLS and bootstrap BCa 95% CIs):\n")
show(slope_table %>% mutate(across(c(anxiety_c_value, slope, SE, ols_low, ols_high, bca_low, bca_high), ~ round(.x, 3)),
                            p = signif(p, 3), interaction_p = signif(interaction_p, 3)))
write_excel_csv(slope_table, out_file("3_simple_slopes.csv"))


# ---- 3.5 Interaction plots ---------------------------------------------------
section("3.5 Interaction plots")

for (name in names(models)) {
  p <- interactions::interact_plot(models[[name]]$step3, pred = presence_c, modx = anxiety_c,
                                   modx.values = "mean-plus-minus", centered = "none",
                                   plot.points = TRUE, point.alpha = 0.4, interval = TRUE,
                                   x.label = "Presence (IPQ mean, centred)",
                                   y.label = outcome_label[[name]],
                                   legend.main = "Anxiety (centred)") +
    ggtitle(paste0(name, ": presence x anxiety")) +
    theme_minimal()
  ggsave(out_file(paste0("3_interaction_", tolower(name), ".png")), p, width = 7, height = 5, dpi = 300)
}
cat("Saved: 3_interaction_rq1.png, 3_interaction_rq2.png\n")


# ---- 3.6 Robustness checks and sensitivity reruns ----------------------------
section("3.6 Robustness checks and sensitivity reruns (final model)")

# Each check refits step 2 and step 3 on its own sample: complete cases on the
# model variables, decision exclusions applied, the listed IDs left out, and
# presence and anxiety re-centred within that sample.
checks <- tribble(
  ~rq,   ~check,    ~label,                                      ~memory,       ~extra,      ~drop,
  "RQ1", "main",    "Main analysis",                             "ravlt_total", "",          NA,
  "RQ1", "A",       "RAVLT delay instead of RAVLT total",        "ravlt_delay", "",          NA,
  "RQ1", "B",       "+ age + sex",                               "ravlt_total", "age + sex", NA,
  "RQ1", "S1",      "without protocol deviation set",            "ravlt_total", "",          "protocol",
  "RQ1", "S2",      "without comprehension / language set",      "ravlt_total", "",          "comprehension",
  "RQ1", "S3",      "without participant characteristics set",   "ravlt_total", "",          "characteristics",
  "RQ1", "A_valid", "RAVLT delay, without RAVLTDelay-only set",  "ravlt_delay", "",          "ravlt_delay_only",
  "RQ2", "main",    "Main analysis",                             "ravlt_total", "",          NA,
  "RQ2", "A",       "RAVLT delay instead of RAVLT total",        "ravlt_delay", "",          NA,
  "RQ2", "B",       "+ age + sex",                               "ravlt_total", "age + sex", NA,
  "RQ2", "C",       "without |z| > 3 (AA01PE, YI16TV)",          "ravlt_total", "",          "rq2_outliers",
  "RQ2", "S1",      "without protocol deviation set",            "ravlt_total", "",          "protocol",
  "RQ2", "S2",      "without comprehension / language set",      "ravlt_total", "",          "comprehension",
  "RQ2", "S3",      "without participant characteristics set",   "ravlt_total", "",          "characteristics",
  "RQ2", "S4",      "without LA19RW (session marked not valid)", "ravlt_total", "",          "la19rw_invalid_session",
  "RQ2", "A_valid", "RAVLT delay, without RAVLTDelay-only set",  "ravlt_delay", "",          "ravlt_delay_only"
)

drop_ids <- function(drop) {
  if (is.na(drop)) character()
  else if (drop == "rq2_outliers") rq2_outlier_ids
  else sensitivity_sets$id[sensitivity_sets$set == drop]
}
rq_outcome <- c(RQ1 = "kol_free_r1", RQ2 = "ssq_change")
rq_excl    <- c(RQ1 = "excl_rq1",    RQ2 = "excl_rq2")

run_check <- function(rq, check, label, memory, extra, drop) {
  extra_vars <- if (extra == "") character() else str_split_1(extra, " \\+ ")
  vars <- c("ipq_mean", "dass_anx", "group", memory, extra_vars, rq_outcome[[rq]])
  s <- d %>%
    filter(if_all(all_of(vars), ~ !is.na(.x)), is.na(.data[[rq_excl[[rq]]]]), !id %in% drop_ids(drop)) %>%
    mutate(y = .data[[rq_outcome[[rq]]]], memory = .data[[memory]],
           presence_c = ipq_mean - mean(ipq_mean), anxiety_c = dass_anx - mean(dass_anx))
  covs <- paste(c("memory", "group", extra_vars), collapse = " + ")
  f2 <- as.formula(paste("y ~", covs, "+ presence_c + anxiety_c"))
  f3 <- as.formula(paste("y ~", covs, "+ presence_c + anxiety_c + presence_c:anxiety_c"))
  m2 <- lm(f2, data = s)
  m3 <- lm(f3, data = s)
  comp <- anova(m2, m3)
  r2_2 <- summary(m2)$r.squared
  r2_3 <- summary(m3)$r.squared

  set.seed(seed)                       # reset before every check
  b <- boot::boot(s, function(dat, idx) coef(lm(f3, data = dat[idx, ])), R = n_boot)

  bca <- map(seq_along(coef(m3)), ~ bca_ci(b, .x))
  coefs <- broom::tidy(m3, conf.int = TRUE) %>%
    mutate(bca_low  = map_dbl(bca, 1),
           bca_high = map_dbl(bca, 2)) %>%
    transmute(rq, check, term, b = estimate, SE = std.error, t = statistic, p = p.value,
              ols_low = conf.low, ols_high = conf.high, bca_low, bca_high)
  inter <- coefs %>% filter(term == "presence_c:anxiety_c")

  cat("\n---", rq, "|", check, "|", label, "| N =", nrow(s),
      "| bootstrap resamples with a missing estimate:", sum(!complete.cases(b$t)), "---\n")
  show(coefs %>% select(-rq, -check) %>%
         mutate(across(c(b, SE, t, ols_low, ols_high, bca_low, bca_high), ~ round(.x, 3)), p = signif(p, 3)))

  list(summary = tibble(rq, check, label, N = nrow(s),
                        n_low = sum(s$group == "low"), n_high = sum(s$group == "high"),
                        n_moderate = sum(s$group == "moderate"),
                        interaction_b = inter$b, SE = inter$SE, p = inter$p,
                        bca_low = inter$bca_low, bca_high = inter$bca_high,
                        delta_R2 = r2_3 - r2_2, F_change = comp$F[2],
                        df1 = comp$Df[2], df2 = comp$Res.Df[2], f2 = (r2_3 - r2_2) / (1 - r2_3),
                        significant_p = inter$p < alpha,
                        bca_excludes_0 = inter$bca_low > 0 | inter$bca_high < 0,
                        n_resamples = b$R),
       coefs = coefs)
}

check_results <- pmap(checks, run_check)
check_summary <- map_dfr(check_results, "summary") %>%
  group_by(rq) %>%
  mutate(same_conclusion_as_main = significant_p == significant_p[check == "main"]) %>%
  ungroup() %>%
  mutate(bca_agrees_with_p = significant_p == bca_excludes_0)

cat("\nSummary: interaction presence x anxiety per check\n")
show(check_summary %>%
       select(rq, check, label, N, interaction_b, p, bca_low, bca_high, delta_R2, F_change, df1, df2, f2,
              same_conclusion_as_main, bca_agrees_with_p) %>%
       mutate(across(c(interaction_b, bca_low, bca_high, F_change, delta_R2, f2), ~ round(.x, 3)),
              p = signif(p, 3)))
cat("same_conclusion_as_main: judged by OLS p. bca_agrees_with_p = FALSE: bootstrap disagrees (report).\n")
write_excel_csv(check_summary, out_file("3_checks_summary.csv"))
check_coefs <- map_dfr(check_results, "coefs")
write_excel_csv(check_coefs, out_file("3_checks_coefficients.csv"))

# Every term in every model where OLS p and the BCa interval disagree (reported
# as borderline; OLS p decides)
mismatches <- check_coefs %>%
  mutate(significant_p = p < alpha, bca_excludes_0 = bca_low > 0 | bca_high < 0) %>%
  filter(significant_p != bca_excludes_0) %>%
  mutate(direction = if_else(significant_p, "OLS significant, BCa includes 0",
                                            "OLS not significant, BCa excludes 0"))
cat("\nOLS p vs BCa disagreements (all terms, all checks):", nrow(mismatches), "\n")
show(mismatches %>% select(rq, check, term, b, p, bca_low, bca_high, direction) %>%
       mutate(across(c(b, bca_low, bca_high), ~ round(.x, 3)), p = signif(p, 3)))
write_excel_csv(mismatches, out_file("3_mismatches.csv"))


# =============================================================================
# PART 4. SELF-CHECK AND SESSION INFO
# =============================================================================

# ---- 4.1 Self-check against the step-by-step analysis ------------------------
section("4.1 Self-check against the logged results")

# Values from the step-by-step analysis (ANALYSIS_LOG Steps 2-7). The script
# stops if a result differs (counts exact; b and p within rounding, .005).
problems <- character()
expect <- function(what, got, want, tol = 0) {
  if (length(got) != length(want) || any(abs(got - want) > tol)) {
    problems <<- c(problems, paste0(what, ": got ", paste(round(got, 4), collapse = " / "),
                                    ", expected ", paste(want, collapse = " / ")))
  }
}

expect("Consenters", nrow(d), 107)
expect("RQ1 / RQ2 / both N", c(nrow(rq1), nrow(rq2), nrow(both)), c(97, 100, 95))
expect("Sex female / male / missing (consenters)",
       c(sum(d$sex == "Female", na.rm = TRUE), sum(d$sex == "Male", na.rm = TRUE), sum(is.na(d$sex))),
       c(74, 32, 1))
bg_all <- background_table %>% select(measure, category, n = `All consenters`)
bg <- function(m, cat) bg_all$n[bg_all$measure == m & bg_all$category == cat]
expect("Currently in education yes / no / missing",
       c(bg("currently in education", "Yes"), bg("currently in education", "No"),
         bg("currently in education", "missing")), c(105, 1, 1))
expect("Nationality Dutch / other / missing / distinct other",
       c(bg("nationality", "Dutch"), bg("nationality", "other"), bg("nationality", "missing"),
         bg("nationality", "distinct other nationalities")), c(56, 50, 1, 26))
expect("Primary language Dutch named / Dutch only",
       c(bg("primary language", "Dutch (alone or with another)"), bg("primary language", "Dutch only")),
       c(56, 54))
expect("FMS peak x SSQ change: n, Pearson r, Spearman rho",
       c(fms_table$n, fms_table$pearson_r, fms_table$spearman_rho), c(99, .751, .716), tol = .0005)
expect("Minimum detectable f2 (RQ1, RQ2)", power_table$min_f2, c(.0827, .0801), tol = .0005)
expect("Cronbach's alpha", alpha_table$alpha, c(.797, .694, .656, .897), tol = .0005)
expect("Sensitivity set sizes",
       as.integer(table(factor(sensitivity_sets$set, levels = c("protocol", "comprehension", "characteristics",
                                                               "ravlt_delay_only", "la19rw_invalid_session")))),
       c(23, 6, 9, 3, 1))
main_inter <- coef_table %>% filter(term == "presence_c:anxiety_c")
expect("Main interaction b (RQ1, RQ2)", main_inter$b, c(-0.200, 0.804), tol = .005)
expect("Main interaction p (RQ1, RQ2)", main_inter$p, c(.0702, .197), tol = .005)
expect("Check N", check_summary$N,
       c(97, 97, 97, 76, 91, 91, 94,  100, 100, 100, 98, 79, 94, 93, 99, 97))
expect("Check interaction b", check_summary$interaction_b,
       c(-0.200, -0.185, -0.198, -0.204, -0.189, -0.256, -0.191,
          0.804,  0.826,  0.810,  0.764,  0.198,  0.812,  1.03,  0.683, 0.814), tol = .005)
expect("Check interaction p", check_summary$p,
       c(.0702, .0976, .071, .091, .0831, .0278, .0951,
         .197, .185, .19, .169, .759, .199, .135, .253, .211), tol = .005)
expect("Checks significant by OLS p (only RQ1 S3)", which(check_summary$significant_p), 6)
expect("Bootstrap resamples (main models, all checks)",
       c(map_int(boots, "R"), check_summary$n_resamples), rep(n_boot, 2 + nrow(check_summary)))

if (length(problems) > 0) {
  cat("SELF-CHECK FAILED:\n", paste("-", problems, collapse = "\n"), "\n")
  stop_here("Results differ from the step-by-step analysis - see above.")
}
cat("Self-check passed: samples, reliability, sensitivity sets, main results and all",
    nrow(check_summary), "checks match the step-by-step analysis.\n")


# ---- 4.2 Session info --------------------------------------------------------
section("4.2 Session info")
print(sessionInfo())
writeLines(capture.output(sessionInfo()), out_file("4_session_info.txt"))
cat("\nAll output is in", out_dir, "\n")

sink()
cat("Results saved to", out_txt, "\n")
