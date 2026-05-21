# Source Code Registry

Use these sources before writing new QC/preprocessing code. Always report the reused source in the final answer.

## Human Endometrium PCOS Atlas

- Paper: `Single-cell profiling of the human endometrium in polycystic ovary syndrome`
- DOI: `10.1038/s41591-025-03592-z`
- Code repository: `https://github.com/ReproductiveEndocrinologyMetabolism/Endo.R`

Reusable files:

- `Spatial RNA-seq analysis/00.1_Stereopy_binning.py`: StereoPy GEF loading, bin-size selection, QC calculation, bin/cell filtering, normalization, clustering checkpoints, h5ad export.

Reusable success points:

- Calculate QC before and after filtering.
- Save a filtered binned h5ad before downstream normalization decisions.
- Spatial QC maps should be checked before interpreting domains or labels.

## GF/SPF Cecum Atlas

- Paper: `Effects of flora deficiency on the structure and function of the large intestine`
- DOI: `10.1016/j.isci.2024.108941`
- Code repository: `https://github.com/1014723815/GF_SPF_cecum`

Reusable files:

- `Spatial transcriptome.R`: GEM-to-Seurat sparse matrix construction, pseudo-image coordinates, BayesSpace and marker spatial plots.

Reusable success points:

- Convert GEM rows into gene x spatial-unit sparse matrices using structured column handling.
- Keep x/y coordinates in metadata and normalize coordinates only for plotting/image compatibility.
- Create coordinate QC outputs before higher-level clustering.

## Cross-Tissue Germ-Free Mouse Atlas

- Paper: `Cross-tissue multi-omics analyses reveal the gut microbiota's absence impacts organ morphology, immune homeostasis, bile acid and lipid metabolism`
- DOI: `10.1002/imt2.272`
- Code repository: `https://github.com/BGI-Intestines/Germ-free-mice`
- Related pipeline evidence: SAW repositories appear repeatedly in local Stereo-seq corpus records.

Reusable files:

- `Figure*/Figure*.ipynb`: figure-level QC and condition plotting after preprocessing.
- `Spleen_*_local.gem.gz`: example GEM-like inputs used by downstream figure notebooks.

Reusable success points:

- Keep per-sample QC tables aligned to sample/condition labels.
- Use consistent coordinate orientation before comparing conditions.
