#!/usr/bin/env python3
"""Describe the archived synthetic event contrast without predicting function."""
import csv
import json
from pathlib import Path
import sys

project = Path(__file__).resolve().parent
root = Path(sys.argv[1]) if len(sys.argv) > 1 else project / "tyk2_focused_rheumatology_immunology/source_results/molecular_qtl_followup/event_annotation"
with (root / "TYK2_grp2_Contained_Structure_Summary.csv").open() as stream:
    structures = list(csv.DictReader(stream))
with (root / "TYK2_Ensembl105_Transcript_Features.csv").open() as stream:
    transcripts = list(csv.DictReader(stream))
selected_id = "ENSG00000105397.grp_2.contained.ENST00000525976"
comparator_id = "ENSG00000105397.grp_2.contained.ENST00000524462"
primary = next(row for row in structures if row["event_id"] == selected_id)
comparator = next(row for row in structures if row["event_id"] == comparator_id)
intervals = lambda row: {tuple(map(int, pair.split("-"))) for pair in row["intervals_ascending_genomic"].split(";")}
missing = intervals(comparator) - intervals(primary)
extra = intervals(primary) - intervals(comparator)
assert missing == {(10351052, 10351162)} and not extra
assert int(comparator["exonic_length"]) - int(primary["exonic_length"]) == 111
scaffold = "ENST00000525621"
mapped_exon = next(row for row in transcripts if row["transcript_id"] == scaffold and row["feature"] == "exon"
                   and row["start"] == "10351052" and row["end"] == "10351162")
variant_exon = next(row for row in transcripts if row["transcript_id"] == scaffold and row["feature"] == "exon"
                    and int(row["start"]) <= 10352442 <= int(row["end"]))
original = next(row for row in transcripts if row["transcript_id"] == "ENST00000525976" and row["feature"] == "transcript")
result = {
    "selected_event": selected_id, "matched_comparator": comparator_id, "reference_release": "Ensembl105/CAGE25",
    "strand": "-", "synthetic_event_exons": 19, "synthetic_event_exonic_bases": 3074,
    "n_peer_event_models": len(structures), "missing_exon_GRCh38": "19:10351052-10351162",
    "missing_exon_bases": 111, "scaffold_transcript": scaffold,
    "scaffold_exon_number": mapped_exon["exon_number"], "variant_scaffold_exon_number": variant_exon["exon_number"],
    "original_label_exon_count": sum(row["transcript_id"] == "ENST00000525976" and row["feature"] == "exon" for row in transcripts),
    "original_transcript_annotation_tag": original["tag"],
    "structural_interpretation": "Exon-24-lacking internal event configuration relative to the named comparator",
    "measurement_interpretation": "Inverse-normalized relative event usage within a group; not exon-24 PSI or absolute isoform abundance",
    "causal_interpretation": "RNA association may accompany the known P1104A protein effect; mediation and causal splicing are unproven"
}
(root / "TYK2_Event_Structural_Interpretation.json").write_text(json.dumps(result, indent=2) + "\n")
with (root / "TYK2_Event_Structural_Interpretation.csv").open("w", newline="") as stream:
    writer = csv.DictWriter(stream, fieldnames=list(result))
    writer.writeheader()
    writer.writerow(result)
print(json.dumps(result))
