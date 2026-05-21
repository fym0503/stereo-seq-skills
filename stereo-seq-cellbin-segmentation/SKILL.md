---
name: stereo-seq-cellbin-segmentation
description: Use when Stereo-seq/STOmics or subcellular spatial transcriptomics data needs cellbin generation, cell segmentation, nuclei/cell masks, ssDNA/DAPI/histology-image alignment, boundary-based expression aggregation, bin-to-cell conversion, segmentation QC, STCellbin, BIDCell, CellSPA, UCS, bin2cell, Thor, or cell-level histology integration.
---

# Stereo-seq Cellbin Segmentation

## Use This For

- Creating or auditing cell-level Stereo-seq expression objects from bin/DNB expression and image-derived cell boundaries.
- Working with cell masks, nuclei masks, ssDNA/DAPI images, tissue masks, image alignment, segmentation QC, or cellbin-to-Seurat/h5ad conversion.
- Comparing or validating segmentation outputs before downstream cell-type mapping, domain discovery, CCI, or histology-linked interpretation.

For expression-only QC after a cellbin object already exists, use `stereo-seq-quality-control-preprocessing`. For downstream label transfer or domains, use the corresponding analysis skill after this step.

## Default Requirements

- Read [source_code.md](references/source_code.md) before designing a workflow. These entries come from real Stereo-seq or subcellular spatial transcriptomics papers and tools with public code.
- If no curated entry fits, search [code_candidates.tsv](references/code_candidates.tsv) for additional article-linked segmentation/histology repositories and reusable files before external search.
- Do not invent cell boundaries from expression alone. Require at least one segmentation source: mask image, boundary polygons, nuclei/cell labels, STCellbin/BIDCell/UCS output, or a clear user-approved segmentation tool.
- Before running Python/R, inspect local environments. Use `conda run -n stereo-skills-py python ...` or `conda run -n stereo-skills-r Rscript ...` when available. Heavy segmentation tools such as STCellbin, BIDCell, Cellpose, Thor, bin2cell, or UCS may require separate installs; if missing, stop that step and tell the user the exact missing package/tool and blocked analysis.
- Keep QC figures paper-ready: Arial, readable labels, equal-aspect spatial/image coordinates, legends outside data, and PDF plus optional 300 dpi PNG.
- In the final response, state the reused paper, DOI, code repository/source file, and dataset-specific edits.

## Workflow

1. Identify inputs:
   - Expression source: GEM/GEF/h5ad/RDS/bin table.
   - Spatial unit: DNB, bin, cellbin, polygon, or image pixel.
   - Segmentation source: cell mask, nuclei mask, boundary polygons, histology/ssDNA/DAPI image, or existing segmentation result.
   - Coordinate transform between expression and image spaces.
2. Read [source_code.md](references/source_code.md) and choose a route by available inputs and closest paper-code evidence, not by hard-coded tissue rules.
3. Use `scripts/cellbin_mask_qc_template.py` for a lightweight mask/coordinate QC report before downstream aggregation or when auditing existing cellbin results.
4. If running a full segmentation method, keep the original tool command/config, model/checkpoint, image channel, pixel size, and post-processing thresholds auditable.
5. Validate segmentation:
   - mask area/cell-size distribution;
   - transcript/bin counts per segmented cell;
   - nuclei/cell boundary overlap if both exist;
   - expression/image coordinate alignment;
   - outlier cells and empty masks.
6. Export cell-level expression/metadata and QC figures before handing off to mapping/domain/CCI skills.

## Reusable Article Code

- `scripts/cellbin_mask_qc_template.py`: lightweight QC and figure template derived from public STCellbin/BIDCell/ascidian cell-segmentation patterns for segmentation summary, mask overlay, and coordinate scatter checks.

When using this or any external method-specific command, tell the user which paper and original source file it came from; use [source_code.md](references/source_code.md) for DOI and repository details.

## Output Expectations

- Segmentation source and coordinate assumptions.
- Cell/mask count, area distribution, transcript/bin count distribution when available.
- QC figures showing mask or boundary geometry and optional expression coordinates.
- Exported metadata table keyed by cell/mask id.
- Blockers for missing images, masks, transforms, packages, or model weights.
- Reused article code source, including paper DOI, code repository/DOI, original file name, and dataset-specific edits.
