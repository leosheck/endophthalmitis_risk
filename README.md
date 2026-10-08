# Endophthalmitis risk

Interactive page estimating the risk of endophthalmitis after intravitreal injection and after cataract surgery, and the factors that change it (prefilled syringe against vial, intracameral antibiotic, and others), with 95% confidence intervals and a numbered source for every figure.

- `index.html`: the page, self-contained (data embedded).
- `data.json`: the extracted series, factor comparisons, pooled estimates and references used by the page.
- `search-log.md`: library and PubMed search log, exclusions and gaps (searched 8 October 2026).
- `src/`: build scripts. `analysis.py` computes exact binomial intervals and DerSimonian-Laird random-effects pools from `ivi_raw.json` and `cat.json` and writes `data.json`; `build_page.py` fills `template.html`. The scripts carry the file paths of the environment they were written in and need those paths edited before rerunning.
