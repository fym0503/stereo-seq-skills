---
name: stereo-seq-publication-story
description: Use whenever a Stereo-seq/STOmics/spatial transcriptomics request asks to reuse, cite, report, compare, or learn from real published Stereo-seq papers, paper templates, article-derived workflows, manuscript figure logic, publication story, scientific story, paper experience, provenance, DOI, code/data source, or "which paper/template was reused". Also use for Chinese requests mentioning 真实论文, 文章模板, 论文模板, story templates, 论文故事线, 复用论文经验, 复用了哪些论文, 原始论文 DOI, 代码/数据来源, 发育图谱, 器官发生, mouse embryo, organogenesis, developmental atlas, or when designing a paper-style analysis plan from a new Stereo-seq dataset. This skill searches bundled unified paper profiles before mapping to local Stereo-seq analysis/plotting skills.
---

# Stereo-seq Publication Story

Use this skill to turn a new Stereo-seq analysis request into a paper-like scientific storyline grounded in previously published Stereo-seq article patterns.

## Use This For

- Rapidly proposing testable scientific questions from a new dataset.
- Matching tissue, species, disease, developmental stage, perturbation, or analysis goal to similar Stereo-seq papers.
- Building a figure-by-figure mini manuscript plan before running downstream workflows.
- Deciding which existing local Stereo-seq skills should be combined, based on similar published analyses rather than fixed tool-selection rules.
- Explaining which paper story patterns were reused and what must still be validated on the user's data.

## Default Requirements

- Search local references first. Do not browse for papers or code unless the bundled corpus is insufficient or the user explicitly asks.
- Start with `scripts/search_paper_stories.py` using the user's tissue, species, condition, stage, biological question, analysis keywords, tool hints, or skill tags. The default search covers all bundled unified paper profiles.
- When code provenance or reusable scripts are needed, also run `scripts/search_code_repositories.py` against the bundled article/method code registry before considering external search.
- Read only the top relevant `references/paper_profiles/papers/Sxxxx.md` files needed for the task.
- Treat the bundled material as a paper-derived evidence base, not as proof about the user's dataset. Separate:
  - what the reference paper showed;
  - what the user's dataset has already demonstrated;
  - what remains a hypothesis to test.
- Do not hard-code tool selection. Infer candidate workflows from similar papers, observed tools, story archetypes, and available local Stereo-seq skills.
- When recommending downstream work, name the relevant local skills and the paper evidence that motivated each choice.
- When writing a final report, include paper id, title, DOI, original code/data source if present in the digest, and whether a story template or digest was reused.

## Workflow

1. Parse the user's dataset and goal: species, tissue, stage, condition, spatial unit, annotations, matched sc/snRNA data, and expected biological question.
2. Search the local paper-story corpus:
   - `python scripts/search_paper_stories.py --query "mouse brain cortex layer spatial domain cell type" --top 8`
   - `python scripts/search_paper_stories.py --species mouse --tissue brain --skill-tags spatial_domain cell_type_mapping spatial_programs --top 8`
3. Search local code-source evidence when executable or figure-template provenance is needed:
   - `python scripts/search_code_repositories.py --query "cellbin segmentation histology boundary Stereo-seq" --top 12`
   - `python scripts/search_code_repositories.py --skill stereo-seq-spatial-domain-discovery --query "mouse embryo domain marker" --top 12`
4. Read the most relevant unified paper profiles first. Use `references/github_code_registry.tsv` and per-skill `references/code_candidates.tsv` only as local source-code navigation, not as a claim that every file was executed.
5. Extract reusable story components:
   - central question and biological gap;
   - figure order template;
   - minimum evidence chain;
   - analysis modules used by the reference paper;
   - validation requirements and common failure modes.
6. Convert the pattern into a dataset-specific plan:
   - candidate scientific questions;
   - analysis modules to run;
   - figure storyboard;
   - expected outputs and decision points;
   - risks, missing metadata, and unsupported claims.
7. Map planned analysis modules to available local skills. Common pairings include:
   - QC/setup: `stereo-seq-quality-control-preprocessing`
   - cellbin/segmentation/histology: `stereo-seq-cellbin-segmentation`
   - annotation/label transfer: `stereo-seq-cell-type-mapping`
   - domain/layer discovery: `stereo-seq-spatial-domain-discovery`
   - marker/pathway/program analysis: `stereo-seq-spatial-programs`
   - spatial interpretation: `stereo-seq-spatial-cell-type-interpretation`
   - CCI: `stereo-seq-spatial-communication` or `stereo-seq-cell-cell-interaction`
   - GRN/regulon: `stereo-seq-spatial-grn-regulon`
   - trajectory/time/state transition: `stereo-seq-developmental-trajectory`
   - 3D/serial sections: `stereo-seq-3d-reconstruction`
   - manuscript figures: `stereo-seq-publication-plotting`
8. Report the provenance of the reused story patterns and code-source templates, then explain how they were adapted to the user's data.

## Reference Layout

- `references/paper_profile_all_index.tsv`: searchable index of all bundled unified paper profiles with tissue, tools, skill tags, story/reuse fields, and profile paths.
- `references/paper_profiles/papers/Sxxxx.md`: one unified profile per paper, combining objective corpus evidence and story/reuse cues.
- `references/story_template_schema.md`: schema used to interpret story templates.
- `references/github_code_registry.tsv`: curated local index of article-owned or method-with-example GitHub/code repositories linked to paper ids, DOI, reusable files, and relevant local skills.
- `scripts/search_code_repositories.py`: searchable interface over the code registry.

## Output Expectations

For a story-design or analysis-planning request, produce:

- 2-5 candidate scientific questions that match the user's dataset.
- A recommended primary story archetype and 1-3 fallback archetypes.
- A concise figure storyboard with required analysis modules.
- Matching reference papers with paper id, title, DOI, and what was reused.
- A mapping from each module to existing local Stereo-seq skills.
- Caveats describing missing metadata, weak evidence, unsupported claims, or package/environment blockers.

For a completed analysis report, include:

- The analyses performed and files generated.
- Which paper-derived profiles were reused.
- Original paper DOI and code/data source as available in the digest.
- What was changed to adapt the template to the user's dataset.
