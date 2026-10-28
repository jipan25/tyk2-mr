# Data sources and reuse

Checked2026-10-28. This repository separates original analytical code from upstream data rights. It includes analysis-derived result tables and limited harmonized statistical values supporting those results, not genome-wide source distributions or participant-level records. No blanket MIT license is asserted over the data directories.

| Source | Attribution and access | Handling in this release |
|---|---|---|
| FinnGen R13 | Kurki et al., Nature2023, doi:10.1038/s41586-022-05473-8; https://www.finngen.fi/en/access_results | Derived MR/colocalization results and selected harmonized values; full regional/genome-wide GWAS and variant-level BF files not mirrored. Follow the official access process for originals. |
| eQTLGen | Vosa et al., Nature Genetics2021, doi:10.1038/s41588-021-00913-z; https://www.eqtlgen.org/cis-eqtls.html | Selected instrument values in analytical tables; full cis/allele-frequency distributions not mirrored. |
| eQTL Catalogue / GTEx / BLUEPRINT | Kerimov et al., Nature Genetics2021, doi:10.1038/s41588-021-00924-w; GTEx doi:10.1126/science.aaz1776; BLUEPRINT doi:10.1016/j.cell.2016.10.026 | Targeted derived molecular results. Catalogue data are CC BY4.0, as explicitly stated at https://www.ebi.ac.uk/eqtl/License/; retain attribution. Catalogue code is Apache2.0 and is not relicensed here. |
| txrevise / Ensembl105 | Alasoo et al., eLife2019, doi:10.7554/eLife.41673; https://doi.org/10.5281/zenodo.6499127; https://doi.org/10.5281/zenodo.7808390 | Extracted event-model/feature summaries, with original source identifiers; no full annotation distribution. |
| OpenGWAS / UKB | https://opengwas.io/; ukb-b-9125; UKB doi:10.1038/s41586-018-0579-z | Targeted summary and harmonization records; no JWT or individual UKB data. Researchers retrieving additional statistics need their own access credentials and must follow upstream terms. |
| 1000 Genomes / SISu4.2 | 1000 Genomes doi:10.1038/nature15393; FinnGen reference documentation | Derived five-variant LD/sensitivity summaries only; no reference-panel genotypes. |
| GSE243917 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE243917; doi:10.1038/s41588-023-01577-7; doi:10.1038/s41588-024-01682-1 | Aggregated gene/cell-compartment summaries; no cell-level matrices or RNA reads. |
| GSE299518 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE299518; doi:10.1038/s41467-025-61842-7 | Derived cell-group detection summaries; no original Seurat object or reads. |

We acknowledge the participants and investigators of FinnGen and the other source studies. Downloadability is not treated as proof of unrestricted redistribution. Users should consult source-specific terms before obtaining or further redistributing originals. The tables are transformations/analytical outputs and must not be described as untouched upstream distributions.

For the72 posterior checks, `coloc_sufficient_statistics.csv` contains three log-sum quantities computed over each shared variant set. Original event selection, allele/variant matching and model assumptions remain documented in the analytical tables. These aggregated inputs cannot reconstruct the full variant-level source datasets. Their input hashes and source identifiers support provenance, not a claim that the complete upstream files are bundled.
