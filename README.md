# TYK2-MR

Code and derived results for a focused analysis of TYK2 molecular signals in rheumatoid arthritis. The related manuscript is **TYK2 molecular signals in rheumatoid arthritis: expression-proxy Mendelian randomization and RNA-event colocalization** (author: Ziqi Yuan; unpublished at this release). The repository account is a hosting/maintenance identity, not a substitute for manuscript authorship. No article DOI or acceptance is claimed.

## Reproduce the core calculations

Use Python3.10-3.12 (tested with Python3.11.14 and NumPy1.26.1; the pinned NumPy release is not intended for newer Python versions):

```sh
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install -r requirements.txt
python3 scripts/reproduce.py
```

The command verifies input hashes and recalculates two FinnGen diagonal IVW estimates, ten leave-one-out estimates, two external-EUR covariance sensitivities, four UKB extraction/harmonization routes, 72 colocalization posterior comparisons and the selected event-model contrast. Outputs are written to `outputs/`; supplied inputs are never modified.

Posterior reproduction uses **derived sufficient statistics**, not redistributed full variant-level FinnGen Bayes factors. It tests the posterior arithmetic but does not repeat variant alignment, fine-mapping fitting or credible-set construction. The original variant-level implementations are retained under `analysis/` for use with upstream inputs; see [ANALYSIS.md](ANALYSIS.md). Do not interpret a successful numerical check as independent biological validation.

## Contents

- `data/tables/`: the28 numbered supplementary tables and Guide, exported as CSV.
- `data/core/`: harmonized analytical values, LD matrix, expected results, event-model summaries and derived Bayes-factor sufficient statistics.
- `data/ukb/`: targeted public-summary extraction/harmonization records for the fixed TYK2 instrument set.
- `data/PROVENANCE.csv`: origin and SHA256 of selected analytical artifacts.
- `data/QTL_SOURCE_URLS.csv`: upstream molecular-QTL locations without local workstation paths.
- `data/EXTERNAL_FINE_MAPPING_INPUTS.json`: hashes of the non-redistributed, target-region fine-mapping inputs used to derive the sufficient statistics.
- `analysis/`: original statistical/annotation implementations, with portable defaults where needed; external inputs and dependencies are required for those broader analyses.

This release contains no manuscript draft, internal review records, author-confirmation workbook, credentials, personal genotypes or raw single-cell matrices. It is not a complete mirror of the original development workspace. Upstream retrieval and original single-cell processing are separate from the offline checks above.

## Scientific interpretation

Expression-proxy estimates are sensitive to LD assumptions and instrument selection. The selected RNA event shares the main RA signal, but signal sharing does not establish RNA mediation or separate RNA effects from P1104A protein effects. Synovial summaries are descriptive: one source pools11RA/oneOA participants and the other analyzed specimen comes from one RA donor. Comparative target counts do not rank efficacy or safety.

## Licensing and attribution

The existing [MIT license](LICENSE) applies to original project code. Third-party data, annotations and statistical values retain their upstream terms; **MIT does not relicense these data**. See [DATA_SOURCES.md](DATA_SOURCES.md) for sources, attribution and redistribution decisions. Upstream files whose redistribution terms were not established are not mirrored as raw datasets.

Codex assisted code development and documentation; human responsibility for scientific interpretation remains with the author. Cite the original FinnGen, eQTLGen, eQTL Catalogue, GTEx, BLUEPRINT and single-cell publications when reusing the scientific results. A release tag identifies the version used for review.
