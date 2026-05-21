# Validation Checks

- Confirm interpretation uses mapped output from a previous method; do not invent labels.
- Check whether enrichment is based on labels, probabilities, or raw counts.
- Use domain-normalized proportions when domains differ in size.
- Use matched bins/cellbins and consistent tissue masks for condition comparisons.
- Preserve uncertainty from deconvolution; avoid hard-label overinterpretation at boundaries.
- Verify coordinate orientation and physical units before making proximity claims.
- Distinguish spatial co-localization from causal interaction.

Interpretation confidence:

- High: replicated across sections/samples and supported by marker/pathology evidence.
- Medium: one section or condition, but consistent with markers and known anatomy.
- Low: no quantitative enrichment, only visual impression.

