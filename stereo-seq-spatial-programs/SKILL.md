---
name: stereo-seq-spatial-programs
description: Use when analyzing Stereo-seq spatial gene programs, including spatially variable genes, spatial co-expression modules, domain or condition DEG, pathway scores, GO/GSEA/GSVA enrichment, TF regulon activity, AUCell/RSS, or NMF/NNMF latent programs.
---

# Stereo-seq Spatial Programs

## Use This For

- Choosing a method for spatially variable genes, spatial modules, DEG, pathway activity, or regulon analysis.
- Interpreting gene programs by domain, condition, time, or cell type in Stereo-seq data.
- Auditing a paper-like program analysis against pilot evidence.

## Default Requirements

- Treat article-derived code reuse as the default when producing analysis code, plots, or tables; the user does not need to ask for this explicitly.
- Before installing packages or creating a new environment, inspect the local Python/R environments. Prefer the existing `stereo-skills-py` and `stereo-skills-r` conda environments when available: `conda run -n stereo-skills-py python ...` or `conda run -n stereo-skills-r Rscript ...`.
- If a required environment or package is missing, stop that step and tell the user exactly what is missing, which command would install it, and which analysis step is blocked. Do not silently install packages, mutate shared environments, or replace the method with an unrelated workaround unless the user asks you to.
- Before writing new plotting or analysis code, read [source_code.md](references/source_code.md) and adapt the closest bundled script in `scripts/`. Do not search GitHub/Zenodo first unless no bundled template fits.
- If no bundled article-derived script fits the task, state that explicitly and keep any custom code minimal.
- In the final response, always state which paper/code source was reused, the paper DOI, the code DOI or repository URL, the original file name, and what was changed for the current dataset.

## Workflow

1. Classify the question:
   - spatial structure of genes,
   - domain/condition DEG,
   - pathway or hallmark scoring,
   - TF/regulon activity,
   - latent factors.
2. Read [tool_selection.md](references/tool_selection.md) as a paper-evidence matrix, then let the question, object stack, desired output, and closest bundled source-code template guide the route.
3. Define outputs before analysis: gene list, module table, DEG table, score matrix, regulon table, enrichment table, or spatial map.
4. Validate with [validation_checks.md](references/validation_checks.md).
5. If code or plotting is needed, read [source_code.md](references/source_code.md) and reuse the matching article-derived script in `scripts/` before writing new code.
6. Ground tool choice in [evidence.md](references/evidence.md).

## Reusable Article Code

- `scripts/p09_hdwgcna_skeleton.R`: adapted from P09 Zenodo `9_Spatial_coexpression.R` for Seurat + hdWGCNA spatial co-expression modules.
- `scripts/gf_cecum_deg_volcano_template.R`: adapted from GF/SPF cecum `Single_cell.R` DEG sections for publication-style volcano plots.
- `scripts/gf_cecum_go_gsea_dotplot_template.R`: adapted from GF/SPF cecum `Single_cell.R` enrichment sections for GO/GSEA dotplots.

When using any bundled script, tell the user which paper and original source file it came from; use [source_code.md](references/source_code.md) for DOI and repository details.

## Output Expectations

- Exact biological question and comparison set.
- Selected tool and why.
- Input object and grouping variables.
- Thresholds and parameters.
- Program output plus spatial/domain interpretation caveats.
- Reused article code source, including paper DOI, code DOI or repository URL, and original file name.
