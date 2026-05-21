# Interpretation Patterns

## Domain or Layer Enrichment

Use when a mapped cell type is enriched in anatomical layers, tissue domains, cancer regions, or disease strata.

Outputs:
- cell type by domain/layer table,
- spatial map with consistent domain colors,
- marker or abundance validation.

## Boundary, Margin, or Interface Localization

Use when the biological claim depends on proximity to an invasive margin, maternal-fetal interface, lesion edge, or tissue boundary.

Outputs:
- distance-to-boundary or layer-stratified summary,
- abundance by boundary layer,
- caveat for coordinate registration and boundary definition.

## Focal Pathology Proximity

Use when interpreting cells around virus-positive spots, stress foci, plaques, tau/A-beta regions, or myocarditic areas.

Outputs:
- radius-band or concentric-level summaries,
- nearest-neighbor/proximity checks,
- control comparison.

## Condition or Time Change

Use when comparing disease/control, treatment, HFD/control, infection time, or developmental stages.

Outputs:
- cell type by condition/time/domain table,
- matched spatial maps with fixed scale,
- explicit note about sample pairing and batch effects.

