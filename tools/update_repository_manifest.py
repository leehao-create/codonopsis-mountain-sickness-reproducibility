#!/usr/bin/env python3
"""Regenerate the repository-wide file manifest and SHA-256 list."""

from __future__ import annotations

import csv
import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "repository_file_manifest.csv"
CHECKSUMS = ROOT / "checksums.sha256"
EXCLUDED = {MANIFEST.resolve(), CHECKSUMS.resolve()}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def classify(relative: str) -> tuple[str, str]:
    parts = relative.split("/")
    module = parts[1] if len(parts) > 2 and parts[0] == "analysis" else "repository"
    if "/scripts/" in relative:
        role = "analysis or rendering script"
    elif "/data/metadata/" in relative:
        role = "sample metadata"
    elif "/data/processed/" in relative:
        role = "processed analysis input"
    elif "/data/dependencies/" in relative:
        role = "frozen analysis dependency"
    elif "/data/derived/" in relative or "/data/reference/" in relative:
        role = "derived or reference input"
    elif "/configs/" in relative:
        role = "analysis configuration"
    elif "/logs/" in relative:
        role = "execution log"
    elif "/results/" in relative:
        role = "frozen final result"
    elif relative.endswith("README.md"):
        role = "documentation"
    else:
        role = "release metadata or module provenance"
    return module, role


def release_status(relative: str) -> tuple[str, str]:
    if relative.startswith("analysis/network_pharmacology/data/derived/Fig1B_") or \
       relative == "analysis/network_pharmacology/data/derived/Fig1C_top15_degree.csv":
        return "yes", "STRING-derived output retained under CC BY 4.0 with attribution in DATA_PROVENANCE.md"
    return "yes", ""


def main() -> None:
    files = sorted(path for path in ROOT.rglob("*") if path.is_file() and path.resolve() not in EXCLUDED)
    rows = []
    checksum_lines = []
    for path in files:
        relative = path.relative_to(ROOT).as_posix()
        digest = sha256(path)
        module, role = classify(relative)
        public_status, note = release_status(relative)
        if module == "WGCNA":
            note = (note + "; " if note else "") + "Authoritative revised beta=30 WGCNA only"
        rows.append({
            "relative_path": relative,
            "bytes": path.stat().st_size,
            "sha256": digest,
            "module": module,
            "role": role,
            "final_or_legacy": "final",
            "public_release": public_status,
            "notes": note,
        })
        checksum_lines.append(f"{digest}  {relative}")

    with MANIFEST.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    CHECKSUMS.write_text("\n".join(checksum_lines) + "\n", encoding="utf-8")
    print(f"Wrote {len(rows)} manifest entries")


if __name__ == "__main__":
    main()
