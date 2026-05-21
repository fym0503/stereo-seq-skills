# Stereo-seq Skills for Codex

Reusable Codex skills for Stereo-seq/STOmics spatial transcriptomics analysis.

This repository is an installable skill bundle. It contains the final skill folders only, not the literature-mining scripts, crawl caches, PDFs, XML files, or intermediate development artifacts used to create them.

## Install

Clone this repository and sync the skill folders into your Codex skills directory:

```bash
git clone https://github.com/fym0503/stereo-seq-skills.git
cd stereo-seq-skills
bash install.sh
```

By default, `install.sh` installs to:

```bash
${CODEX_HOME:-$HOME/.codex}/skills
```

You can override the target:

```bash
CODEX_HOME=/path/to/.codex bash install.sh
```

Restart Codex after installation so the new skills are discovered.

## Included Skills

- `stereo-seq-quality-control-preprocessing`
- `stereo-seq-publication-plotting`
- `stereo-seq-publication-story`
- `stereo-seq-spatial-domain-discovery`
- `stereo-seq-spatial-programs`
- `stereo-seq-spatial-cell-type-interpretation`
- `stereo-seq-cell-type-mapping`
- `stereo-seq-cell-cell-interaction`
- `stereo-seq-spatial-communication`
- `stereo-seq-spatial-grn-regulon`
- `stereo-seq-developmental-trajectory`
- `stereo-seq-3d-reconstruction`

## Publication Story Corpus

The `stereo-seq-publication-story` skill includes bundled local references:

- 462 objective Stereo-seq paper digests.
- 50 reusable paper story templates.
- Search script: `stereo-seq-publication-story/scripts/search_paper_stories.py`.

The skill uses these local references to help Codex identify similar published Stereo-seq analyses, reuse paper-style figure logic, and report DOI/source provenance when adapting workflows to a new dataset.

## Environment

These skills describe workflows and include reusable scripts/templates, but they do not install Python or R packages automatically. If a workflow needs unavailable packages, Codex should report the missing packages and ask the user to install or activate a suitable environment.

Typical analysis environments may need packages such as `scanpy`, `anndata`, `squidpy`, `stereopy`, `GraphST`, `pyscenic`, `CellChat`, `Seurat`, or related method-specific dependencies depending on the requested workflow.
