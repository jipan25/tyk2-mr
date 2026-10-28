# Original analysis implementations

The supported offline entry point is `python3 scripts/reproduce.py`. It uses only the curated inputs in this repository. The scripts below retain the project's original statistical implementations for additional inspection/reanalysis; they are **not** a promise that all original source datasets are present.

| Script | Purpose | Additional inputs |
|---|---|---|
| j63 | SNP Q contributions and leave-one-out influence | Harmonized clinical records and archived MR outputs in the layout read by the script |
| j65 | Complete-region gene ABF, MR and sensitivity estimates | FinnGen R13 regional summary statistics/manifest; full eQTLGen TYK2 cis records and allele frequencies; five-instrument RDS |
| j72 | Signal-specific coloc from published SuSiE BFs | FinnGen regional BFs/credible sets; four Catalogue event-BF and credible-set files; selected metadata |
| j74 / j78 | Annotation extraction and event-structure interpretation | Original txrevise/Ensembl archives (j74 lists exact URLs); extracted feature/model tables |
| j75 | Variant-level independent IVW and posterior checks | The full target-region shared-variant BF inputs, unlike the sufficient-statistic check in the supported entry point |
| j79 | External-EUR covariance sensitivity | Five-variant signed LD and harmonized records |
| j84 / j85 / j86 | GTEx-lead, regional ABF and fixed-marker checks | Original GTEx leads, FinnGen regional GWAS, six Catalogue gene-level regional extracts |

The original R analyses used R4.3.3, TwoSampleMR0.7.4, MRPRESSO1.0 and coloc5.2.3, with data.table, ggplot2 and jsonlite. The Python numerical/annotation scripts primarily use the standard library; j79 also requires NumPy. R scripts accept their input-root arguments at the top of each file. Historical workstation defaults were changed to relative `data/upstream/...` paths; statistical formulas were not changed. Some Python scripts retain the historical project-relative source layout under their script directory; pass explicit paths where supported or reproduce that layout when performing the broader source-level checks.

Obtain FinnGen R13 through https://www.finngen.fi/en/access_results and its linked release documentation. The summary manifest is https://storage.googleapis.com/finngen-public-data-r13/summary_stats/finngen_R13_manifest.tsv. The clinical regional extraction used GRCh38 chr19:9365604-11365604; the fine-mapping signal region was chr19:8852442-11852442. These are different inputs and must not be substituted for one another. Follow the release manifest/documentation for GWAS and SuSiE file URLs.

Catalogue URLs and dataset IDs are in `data/QTL_SOURCE_URLS.csv`. The selected event is `ENSG00000105397.grp_2.contained.ENST00000525976`; the four event contexts are QTD000024, QTD000034, QTD000329 and QTD000359. The six gene contexts are described by the supplementary availability tables. Use the source's chromosome/build conventions and archived metadata when extracting regions; do not silently substitute current releases. Source hashes for the target-region BF extracts are recorded separately; independently recompressed GZIP bytes can differ, so compare scientific contents as well as byte-level provenance.

The repository does not rerun raw single-cell processing, a new SuSiE fit, the historical full OpenGWAS comparator search or manuscript/PI assembly. All28 published-analysis tables are provided to inspect the wider results. Replicating broader steps requires the source datasets and original dependencies listed above.
