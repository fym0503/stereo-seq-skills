#!/usr/bin/env python3
"""Search bundled Stereo-seq article/method code repositories."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path
from typing import Iterable


TOKEN_RE = re.compile(r"[A-Za-z0-9]+")


def tokenize(text: str) -> set[str]:
    return {token.lower() for token in TOKEN_RE.findall(text or "") if len(token) > 1}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def split_terms(values: Iterable[str] | None) -> set[str]:
    terms: set[str] = set()
    for value in values or []:
        for part in re.split(r"[;,]\s*|\s+", value):
            part = part.strip().lower()
            if part:
                terms.add(part)
    return terms


def contains_term(haystack: str, term: str) -> bool:
    return term.lower() in (haystack or "").lower()


def score_record(record: dict[str, str], args: argparse.Namespace) -> float:
    score = 0.0
    weighted_fields = {
        "repo": 4.0,
        "paper_titles": 3.0,
        "paper_dois": 2.0,
        "description": 2.0,
        "skill_names": 3.0,
        "raw_skill_tags": 2.0,
        "reusable_files": 2.0,
        "tissue": 1.5,
        "species": 1.5,
    }
    query_tokens = tokenize(args.query)
    for field, weight in weighted_fields.items():
        score += len(query_tokens & tokenize(record.get(field, ""))) * weight

    try:
        score += min(float(record.get("score", "0")), 900.0) / 100.0
    except ValueError:
        pass

    for term in split_terms(args.skill):
        if term in record.get("skill_names", "").lower():
            score += 8.0
    for term in split_terms(args.paper_id):
        if term in record.get("linked_papers", "").lower():
            score += 20.0
    for term in split_terms(args.repo):
        if term in record.get("repo", "").lower():
            score += 20.0

    query_text = " ".join(sorted(query_tokens))
    if any(term in query_text for term in {"cellbin", "segmentation", "histology", "boundary", "mask"}):
        cellbin_terms = [
            "stereo-seq-cellbin-segmentation",
            "cellbin",
            "segmentation",
            "cellseg",
            "cell mask",
            "nuclei",
            "histology",
            "boundary",
            "bin2cell",
            "bidcell",
            "stcellbin",
            "cellspa",
            "thor",
        ]
        haystack = " ".join(
            [
                record.get("repo", ""),
                record.get("description", ""),
                record.get("skill_names", ""),
                record.get("reusable_files", ""),
            ]
        ).lower()
        score += sum(4.0 for term in cellbin_terms if term in haystack)
    return score


def passes_filters(record: dict[str, str], args: argparse.Namespace) -> bool:
    if args.skill:
        haystack = record.get("skill_names", "").lower()
        if not all(term in haystack for term in split_terms(args.skill)):
            return False
    if args.paper_id:
        haystack = record.get("linked_papers", "").lower()
        if not all(term in haystack for term in split_terms(args.paper_id)):
            return False
    if args.repo:
        haystack = record.get("repo", "").lower()
        if not all(contains_term(haystack, term) for term in split_terms(args.repo)):
            return False
    if args.classification:
        haystack = record.get("classification", "").lower()
        if not all(term in haystack for term in split_terms(args.classification)):
            return False
    return True


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Search local Stereo-seq GitHub/code-source registry."
    )
    parser.add_argument("--query", default="", help="Free-text search query.")
    parser.add_argument("--skill", nargs="*", help="Required skill name or term.")
    parser.add_argument("--paper-id", nargs="*", help="Required paper id, such as S0097.")
    parser.add_argument("--repo", nargs="*", help="Required repository name term.")
    parser.add_argument("--classification", nargs="*", help="Classification filter.")
    parser.add_argument("--top", type=int, default=12, help="Number of records to show.")
    args = parser.parse_args()

    ref_dir = Path(__file__).resolve().parent.parent / "references"
    registry_path = ref_dir / "github_code_registry.tsv"
    if not registry_path.exists():
        raise FileNotFoundError(f"Missing code registry: {registry_path}")

    rows = [row for row in read_tsv(registry_path) if passes_filters(row, args)]
    scored = [(score_record(row, args), row) for row in rows]
    scored.sort(key=lambda item: (-item[0], item[1].get("repo", "").lower()))

    fields = [
        "score",
        "repo",
        "classification",
        "linked_papers",
        "paper_dois",
        "paper_titles",
        "skill_names",
        "source_url",
        "reusable_files",
    ]
    print("\t".join(fields))
    for score, row in scored[: max(0, args.top)]:
        out = {
            "score": f"{score:.2f}",
            "repo": row.get("repo", ""),
            "classification": row.get("classification", ""),
            "linked_papers": row.get("linked_papers", ""),
            "paper_dois": row.get("paper_dois", ""),
            "paper_titles": row.get("paper_titles", ""),
            "skill_names": row.get("skill_names", ""),
            "source_url": row.get("source_url", ""),
            "reusable_files": row.get("reusable_files", ""),
        }
        print("\t".join(out[field].replace("\t", " ") for field in fields))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
