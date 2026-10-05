# Network pharmacology

This module contains compact derived outputs used to render Figure 1 and the final 54-gene candidate list. It does not rerun the historical database searches. Raw GeneCards, OMIM, TCMSP and HERB exports are intentionally excluded.

The three frozen STRING-derived tables are retained under STRING's Creative Commons Attribution 4.0 license. STRING states that this license covers tables, data files, individual scores and data points obtained from its website, APIs or downloads. Source: <https://version-12-5.string-db.org/cgi/access?footer_active_subpage=licensing>. Cite STRING and indicate that the node metrics/plot tables are derived and reformatted.

The archived network contains 54 human nodes and 543 edges with minimum combined score 0.400. Recoverable query parameters are: species `9606`, functional network, required score `400`, and no added interactors. The exact historical STRING release was not recorded, so current database output is not expected to be byte-identical. `scripts/retrieve_current_STRING_network.py` documents a version-specific STRING v12.5 API request and writes current results to `reconstructed_data/` without overwriting frozen Figure 1 tables.

Frozen table checksums are preserved in `data/derived/Fig1_frozen_source_manifest.csv` and the repository manifest. `scripts/draw_Fig1_final_candidate.py` consumes the frozen CC BY 4.0 tables and does not query STRING or recalculate the archived enrichment results.
