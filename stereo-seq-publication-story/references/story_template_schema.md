# Stereo-seq Paper Story Template Schema

Purpose: turn a Stereo-seq paper into a reusable story pattern, not just a list of methods. This schema is intentionally open-ended for method modules; only the fields are fixed.

## Paper Metadata

- `paper_id`
- `title`
- `journal`
- `year`
- `doi`
- `code_urls`
- `data_urls`
- `evidence_source`: e.g. full text, figure captions, methods, code, manual reading
- `review_status`: draft / manually_audited / code_verified

## Story Layer

- `central_question`: What is the biological or technical question?
- `biological_gap`: What was unknown before this paper?
- `technical_gap`: Why was Stereo-seq or spatial multi-omics needed?
- `one_sentence_story`: The paper's story in one sentence.
- `story_archetype`: atlas, anatomical_domain, developmental_timecourse, disease_niche, regeneration, spatial_mechanism, method, cross_dataset_trait_mapping, 3d_4d_resource, validation_translation, custom.
- `final_claim_type`: resource, descriptive atlas, mechanism hypothesis, experimentally validated mechanism, clinical/trait association, method benchmark.
- `main_axis`: anatomy, time, condition, disease state, cell state, spatial boundary, perturbation, genotype, trait, cross-species, custom.

## Dataset Layer

- `species`
- `tissue_or_system`
- `conditions_or_timepoints`
- `n_sections_or_samples`
- `spatial_unit`: DNB, bin, spot, cellbin, segmented cell, pseudo-spot, unknown.
- `matched_modalities`: scRNA/snRNA, scATAC, histology, IF/IHC/FISH, metabolomics, genetics/GWAS, public datasets, perturbation, clinical outcomes.
- `validation_data`

## Figure Logic

Each figure should be summarized as:

- `figure`
- `role`: setup, atlas, spatial_structure, cell_type_mapping, program, interaction, trajectory, regulatory_network, validation, resource, method, translation, custom.
- `claim_supported`
- `input_data`
- `analysis`
- `plot_types`
- `key_result`
- `transition_to_next`
- `reusable_pattern`
- `caveat`

## Open Method Modules

Use one entry for every logical method block. Do not force a fixed list.

- `module_name`
- `module_type`: preprocessing, annotation, spatial_structure, differential_program, interaction, trajectory, regulatory_network, morphology, validation, cross_dataset_comparison, resource_building, statistics, method_modeling, clinical_association, custom.
- `biological_question`
- `why_used_here`
- `input_data`
- `output_artifacts`
- `tools_or_code`
- `figures_supported`
- `claim_supported`
- `evidence_strength`: exploratory, descriptive, statistical, replicated, experimentally_validated.
- `transferability`
  - `reusable_when`
  - `not_reusable_when`
- `implementation_notes`

## Extracted Reusable Story Pattern

- `figure_order_template`
- `minimum_evidence_chain`: e.g. spatial map -> marker DEG -> pathway -> validation.
- `claim_templates`: reusable claim wordings with required evidence.
- `analysis_modules_to_call`: names of existing local Stereo-seq skills to combine.
- `transferable_to_new_data`
- `not_transferable_without`
- `common_failure_modes`

## Audit Notes

- `what_is_inferred_from_captions`
- `what_requires_full_methods_check`
- `what_requires_code_check`
- `what_requires_domain_expert_review`
