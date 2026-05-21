# Source Code Registry

Use these sources when adapting spatial interpretation figures. Always report the reused source in the final answer.

## P09 AD Prefrontal Cortex

- Paper: `Stereo-seq of the prefrontal cortex in aging and Alzheimer's disease`
- DOI: `10.1038/s41467-024-54715-y`
- Code record: Zenodo `10.5281/zenodo.14048103`
- Source record: `https://zenodo.org/records/14048103`

Reusable files:

- `5_location_Tau.py`: target-cell spatial localization with highlighted target cells in red and background cells in grey.
- `case_tau_location.py` and `control_tau_location.py`: condition-stratified target-cell localization after subtype annotation.
- `2_layer_annotation.ipynb`: layer/domain proportion pies and matched spatial maps.

Reusable success points:

- Highlight focal/pathology-associated cells in red against grey background cells.
- Use the same target/background legend across raw and QC-filtered spatial maps.
- Invert coordinate axes deliberately and document the orientation.
- Use per-condition matched panels for case/control interpretation.

## GF/SPF Cecum Atlas

- Paper: `Effects of flora deficiency on the structure and function of the large intestine`
- DOI: `10.1016/j.isci.2024.108941`
- Code repository: `https://github.com/1014723815/GF_SPF_cecum`

Reusable files:

- `Single_cell.R`: condition-wise cell-type proportion summaries and comparison plots.

Bundled template:

- `scripts/gf_cecum_celltype_proportion_template.R`: stacked/fraction barplot of mapped cell types by condition or sample.

Reusable success points:

- Use stable cell-type color ordering across condition panels.
- Keep legends large enough to read in paper figures.
- Export both sample-level and condition-level proportions when sample labels are available.
