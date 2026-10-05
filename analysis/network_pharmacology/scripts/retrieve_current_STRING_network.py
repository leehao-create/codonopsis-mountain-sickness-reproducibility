#!/usr/bin/env python3
"""Retrieve a current STRING network without overwriting the frozen release tables.

This helper documents the recoverable query parameters. STRING is versioned and
updated, so current output is not expected to be byte-identical to the archived
Figure 1 inputs.
"""

from pathlib import Path
import csv
import urllib.parse
import urllib.request

import networkx as nx


MODULE = Path(__file__).resolve().parents[1]
GENES = MODULE / "data" / "derived" / "Codonopsis_mountain_sickness_54_targets.csv"
OUT = MODULE / "reconstructed_data"
API = "https://version-12-5.string-db.org/api/tsv/network"


def main() -> None:
    with GENES.open(newline="", encoding="utf-8-sig") as stream:
        reader = csv.reader(stream)
        next(reader)
        genes = [row[0].strip() for row in reader if row and row[0].strip()]
    if len(genes) != 54 or len(set(genes)) != 54:
        raise ValueError("Expected the frozen 54-gene candidate list")

    payload = urllib.parse.urlencode({
        "identifiers": "\r".join(genes),
        "species": 9606,
        "network_type": "functional",
        "required_score": 400,
        "add_nodes": 0,
        "caller_identity": "Codonopsis_mountain_sickness_reproducibility_repository",
    }).encode("utf-8")
    request = urllib.request.Request(API, data=payload, method="POST")
    with urllib.request.urlopen(request) as response:
        content = response.read().decode("utf-8")

    OUT.mkdir(parents=True, exist_ok=True)
    edge_file = OUT / "STRING_v12.5_54_targets_required_score_400.tsv"
    edge_file.write_text(content, encoding="utf-8")

    rows = list(csv.DictReader(content.splitlines(), delimiter="\t"))
    graph = nx.Graph()
    graph.add_nodes_from(genes)
    graph.add_edges_from((row["preferredName_A"], row["preferredName_B"]) for row in rows)
    degree_file = OUT / "STRING_v12.5_54_targets_node_degree.csv"
    with degree_file.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["gene", "degree"])
        writer.writerows(sorted(graph.degree, key=lambda item: (-item[1], item[0])))

    print(f"Wrote {edge_file}")
    print(f"Wrote {degree_file}")
    print("Current STRING results may differ from the frozen Figure 1 tables.")


if __name__ == "__main__":
    main()
