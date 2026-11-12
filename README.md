# Does Anxiety Moderate Presence–Memory and Presence–Cybersickness Relationships? A Virtual Reality Method of Loci Study

**Author:** Ana Carolina Martins Tramujas\
Bachelor thesis, BSc Cognitive Science and Artificial Intelligence, Tilburg University\
Supervisor: Dr. Eriko Fukuda\
Second reader: Sasha Kanjeeva

## What this script does

`scripts/vrmol_final_analysis.R` runs the full analysis, from the raw data
file to the final results. It reads and cleans the data, computes the scale
scores, runs the descriptives and assumption checks, fits two hierarchical
moderated regressions with bootstrap (BCa) confidence intervals and
robustness checks, and finally checks its own results against the
step-by-step analysis.

- **RQ1:** Does anxiety moderate the relationship between sense of presence
  and memory performance, controlling for baseline memory and immersion group?
- **RQ2:** Does anxiety moderate the relationship between sense of presence
  and cybersickness severity, controlling for baseline memory and immersion group?

## Software

R 4.4.2, with these packages:

| Package | Version |
|---|---|
| readr | 2.2.0 |
| dplyr | 1.2.1 |
| tidyr | 1.3.1 |
| stringr | 1.5.1 |
| purrr | 1.0.4 |
| tibble | 3.2.1 |
| ggplot2 | 4.0.3 |
| car | 3.1-3 |
| lmtest | 0.9-40 |
| interactions | 1.2.0 |
| psych | 2.4.6.26 |
| boot | 1.3-31 |
| broom | 1.0.8 |

The full session information is printed at the end of each run.

## How to run

1. Place the raw data file at `data/raw/kingdom_of_loci_raw.csv`.
2. From the project folder, run:

   ```
   Rscript scripts/vrmol_final_analysis.R
   ```

The results are written to `output/final/`. A run takes about 5 to 10
minutes (5000 bootstrap resamples, fixed seed 20261012).

## Data availability

The data belong to the VRMoL feasibility study team at Tilburg University
and are not publicly available. This repository contains no data.
