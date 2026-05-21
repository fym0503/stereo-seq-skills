---
name: stereo-seq-spatial-cell-type-interpretation
description: Use after Stereo-seq cell types or states have been mapped to spatial coordinates, when the task is to interpret where cell populations are enriched, depleted, layered, associated with tissue domains, pathology, interfaces, treatments, or time points, and to turn those spatial patterns into cautious biological conclusions.
---

# Stereo-seq Spatial Cell-Type Interpretation

## Use This For

- Interpreting mapped cell-type labels, probabilities, or abundances in tissue coordinates.
- Comparing cell-type localization by domain, layer, pathology, interface, treatment, condition, or time.
- Producing concise biological interpretations with supporting evidence and caveats.

Do not use this skill to choose or run deconvolution; use `stereo-seq-cell-type-mapping` first.

## Default Requirements

- Treat article-derived code reuse as the default when producing analysis code, plots, or tables; the user does not need to ask for this explicitly.
- Before installing packages or creating a new environment, inspect the local Python/R environments. Prefer the existing `stereo-skills-py` and `stereo-skills-r` conda environments when available: `conda run -n stereo-skills-py python ...` or `conda run -n stereo-skills-r Rscript ...`.
- If a required environment or package is missing, stop that step and tell the user exactly what is missing, which command would install it, and which analysis step is blocked. Do not silently install packages, mutate shared environments, or replace the method with an unrelated workaround unless the user asks you to.
- Before writing new plotting or analysis code, read [source_code.md](references/source_code.md) and adapt the closest bundled script in `scripts/`. Do not search GitHub/Zenodo first unless no bundled template fits.
- If no bundled article-derived script fits the task, state that explicitly and keep any custom code minimal.
- In the final response, always state which paper/code source was reused, the paper DOI, the code DOI or repository URL, the original file name, and what was changed for the current dataset.

## Workflow

1. Confirm required inputs:
   - Spatial coordinates.
   - Cell-type labels or abundance/probability matrix.
   - Domain/layer/region/condition/time labels if the question involves comparisons.
   - Optional pathology, histology, infection, stress, margin, or interface annotations.
2. Summarize localization:
   - cell type by domain/layer/region,
   - enrichment/depletion near boundaries or foci,
   - condition/time changes,
   - uncertainty from mapping.
3. Use [interpretation_patterns.md](references/interpretation_patterns.md) to match the biological question.
4. Validate claims with [validation_checks.md](references/validation_checks.md).
5. If code or plotting is needed, read [source_code.md](references/source_code.md) and reuse the matching article-derived script in `scripts/` before writing new code.
6. Ground the explanation in [evidence.md](references/evidence.md) when choosing an interpretation pattern.

## Reusable Article Code

- `scripts/p09_target_cell_location.py`: adapted from P09 Zenodo `5_location_Tau.py`, `case_tau_location.py`, and `control_tau_location.py` for red/grey target-cell spatial localization plots.
- `scripts/gf_cecum_celltype_proportion_template.R`: adapted from GF/SPF cecum `Single_cell.R` for condition-wise cell-type proportion plots after mapping.

When using any bundled script, tell the user which paper and original source file it came from; use [source_code.md](references/source_code.md) for DOI and repository details.

## Output Expectations

Return:

- Main spatial pattern.
- Quantitative support needed or already available.
- Biological interpretation in cautious language.
- Alternative explanations and caveats.
- Figures/tables that should support the claim.
- Reused article code source, including paper DOI, code DOI or repository URL, and original file name.
