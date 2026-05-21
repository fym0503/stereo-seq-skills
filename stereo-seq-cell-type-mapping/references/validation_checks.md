# Validation Checks

## Required Checks

- Confirm gene identifier compatibility between spatial data and reference.
- Check reference tissue, species, disease state, and developmental stage.
- For mixed bins, keep probability or abundance values, not only the max label.
- Validate top labels with marker genes on spatial maps.
- Inspect ambiguous or low-confidence bins near tissue boundaries.
- Compare cell-type calls to existing domains, layers, histology, or known anatomy.

## Warning Signs

- A cell type absent from the reference dominates large regions.
- Marker genes disagree with the assigned cell type.
- All bins receive a single dominant cell type despite known mixed tissue.
- Rare cell types are called in low-UMI or high-background areas.
- Coordinate flips or rotations make marker localization biologically implausible.

## Reporting

Include a short confidence statement:

- High: reference matches tissue/stage, markers agree, spatial pattern is coherent.
- Medium: reference is plausible but incomplete or marker evidence is mixed.
- Low: marker-only labels, mismatched reference, or high ambiguity.

