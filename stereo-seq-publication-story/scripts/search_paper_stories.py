#!/usr/bin/env python3
"""Search bundled Stereo-seq paper-story templates and paper digests."""

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


def merged_records(ref_dir: Path, digest_index: str) -> list[dict[str, str]]:
    digest_path = ref_dir / digest_index
    if not digest_path.exists():
        raise FileNotFoundError(f"Missing digest index: {digest_path}")
    digest_rows = read_tsv(digest_path)
    story_rows = read_tsv(ref_dir / "story_template_50_index.tsv")
    story_by_id = {row["paper_id"]: row for row in story_rows}

    records: list[dict[str, str]] = []
    for row in digest_rows:
        record = dict(row)
        story = story_by_id.get(row["paper_id"], {})
        for key, value in story.items():
            if key in {"selection_rank", "doi", "year", "journal", "title", "paper_id"}:
                continue
            record[key] = value
        records.append(record)
    return records


def score_record(record: dict[str, str], args: argparse.Namespace) -> float:
    score = 0.0
    weighted_fields = {
        "title": 3.0,
        "tissue": 3.0,
        "species": 2.0,
        "story_archetype": 3.0,
        "main_axis": 3.0,
        "final_claim_type": 2.0,
        "skill_tags": 3.0,
        "tools_observed": 2.0,
        "journal": 0.5,
        "year": 0.5,
    }

    query_tokens = tokenize(args.query)
    for field, weight in weighted_fields.items():
        field_tokens = tokenize(record.get(field, ""))
        score += len(query_tokens & field_tokens) * weight

    if args.paper_id and record.get("paper_id", "").lower() == args.paper_id.lower():
        score += 100.0
    if args.species and contains_term(record.get("species", ""), args.species):
        score += 8.0
    if args.tissue and contains_term(record.get("tissue", ""), args.tissue):
        score += 8.0

    skill_tags = record.get("skill_tags", "").lower()
    tools = record.get("tools_observed", "").lower()
    archetype = record.get("story_archetype", "").lower()

    for term in split_terms(args.skill_tags):
        if term in skill_tags:
            score += 6.0
    for term in split_terms(args.tools):
        if term in tools:
            score += 5.0
    for term in split_terms(args.archetype):
        if term in archetype:
            score += 5.0

    try:
        rank = int(record.get("selection_rank", "999"))
    except ValueError:
        rank = 999
    score += max(0.0, 3.0 - rank / 25.0)
    return score


def passes_filters(record: dict[str, str], args: argparse.Namespace) -> bool:
    if args.paper_id and record.get("paper_id", "").lower() != args.paper_id.lower():
        return False
    if args.species and not contains_term(record.get("species", ""), args.species):
        return False
    if args.tissue and not contains_term(record.get("tissue", ""), args.tissue):
        return False
    if args.skill_tags:
        skill_tags = record.get("skill_tags", "").lower()
        if not all(term in skill_tags for term in split_terms(args.skill_tags)):
            return False
    if args.tools:
        tools = record.get("tools_observed", "").lower()
        if not all(term in tools for term in split_terms(args.tools)):
            return False
    if args.archetype:
        archetype = record.get("story_archetype", "").lower()
        if not all(term in archetype for term in split_terms(args.archetype)):
            return False
    return True


def relative_path(record: dict[str, str], key: str) -> str:
    value = record.get(key, "")
    if key == "digest_path" and value:
        return "references/" + value
    if key == "template_path" and value:
        return "references/" + value
    return value


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Search local Stereo-seq publication story templates and paper digests."
    )
    parser.add_argument("--query", default="", help="Free-text search query.")
    parser.add_argument("--paper-id", help="Exact paper id, such as S0097.")
    parser.add_argument("--species", help="Species filter, such as mouse or human.")
    parser.add_argument("--tissue", help="Tissue/system filter, such as brain or placenta.")
    parser.add_argument("--skill-tags", nargs="*", help="Required skill tag terms.")
    parser.add_argument("--tools", nargs="*", help="Required observed tool terms.")
    parser.add_argument("--archetype", nargs="*", help="Required story archetype terms.")
    parser.add_argument("--top", type=int, default=8, help="Number of records to show.")
    parser.add_argument(
        "--digest-index",
        default="paper_digest_all_index.tsv",
        help="Digest index in references/. Defaults to the all-paper index.",
    )
    args = parser.parse_args()

    script_dir = Path(__file__).resolve().parent
    ref_dir = script_dir.parent / "references"
    records = merged_records(ref_dir, args.digest_index)

    scored = [
        (score_record(record, args), record)
        for record in records
        if passes_filters(record, args)
    ]
    scored.sort(key=lambda item: (-item[0], int(item[1].get("selection_rank", "999"))))

    fields = [
        "score",
        "paper_id",
        "doi",
        "year",
        "journal",
        "title",
        "species",
        "tissue",
        "story_archetype",
        "main_axis",
        "final_claim_type",
        "skill_tags",
        "tools_observed",
        "story_template",
        "paper_digest",
    ]
    # Use print instead of pandas so the script works in a minimal Python env.
    print("\t".join(fields))
    for score, record in scored[: max(0, args.top)]:
        row = {
            "score": f"{score:.2f}",
            "paper_id": record.get("paper_id", ""),
            "doi": record.get("doi", ""),
            "year": record.get("year", ""),
            "journal": record.get("journal", ""),
            "title": record.get("title", ""),
            "species": record.get("species", ""),
            "tissue": record.get("tissue", ""),
            "story_archetype": record.get("story_archetype", ""),
            "main_axis": record.get("main_axis", ""),
            "final_claim_type": record.get("final_claim_type", ""),
            "skill_tags": record.get("skill_tags", ""),
            "tools_observed": record.get("tools_observed", ""),
            "story_template": relative_path(record, "template_path"),
            "paper_digest": relative_path(record, "digest_path"),
        }
        print("\t".join(row[field].replace("\t", " ") for field in fields))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
