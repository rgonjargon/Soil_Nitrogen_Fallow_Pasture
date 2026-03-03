# Fate of soil nitrogen under long-term bare fallow

Analysis comparing total soil N, POM, microbial biomass, and cations between permanent pasture and permanent fallow over time. Fitted with mixed-effects models (lme4) and model selection via AICc.

## Project structure

| Folder / file | Contents |
|---------------|----------|
| **Scripts/** | `1_analysis.qmd` — main Quarto report; renders to `Analysis_Report.html` |
| **Analysis/new_data/** | Input data (Excel) used by the report |

## How to run

From the project root (`230601_Trish_Fallow_vs_Pasture`):

```bash
quarto render Scripts/1_analysis.qmd
```

Output: `Scripts/Analysis_Report.html`.

## Requirements

- [R](https://www.r-project.org/) and [Quarto](https://quarto.org/)
- R packages: see the setup chunk in `1_analysis.qmd` (tidyverse, readxl, here, lme4, MuMIn, performance, janitor, RColorBrewer, gridExtra, colorspace, etc.)
- CmdStan installed and path set in the setup chunk if using brms

## Author

Dr Tom Moore — Statistical Scientist
