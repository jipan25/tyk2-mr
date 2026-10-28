#!/usr/bin/env python3
"""Extract frozen event and transcript structures; do not infer RNA mediation."""
import argparse
import csv
import gzip
import hashlib
import json
import re
import tarfile
import time
import urllib.request
from pathlib import Path

GENE = "ENSG00000105397"
EVENT = GENE + ".grp_2.contained.ENST00000525976"
VARIANT_POSITION = 10352442
SOURCES = {
    "Homo_sapiens.GRCh38.105.txrevise.CAGE25.tar.gz": (
        "https://zenodo.org/records/6499127/files/Homo_sapiens.GRCh38.105.txrevise.CAGE25.tar.gz?download=1",
        "17add01f40f1584f675d58a15a16baec"),
    "txrevise_Ensembl_105_phenotype_metadata.tsv.gz": (
        "https://zenodo.org/records/7808390/files/txrevise_Ensembl_105_phenotype_metadata.tsv.gz?download=1",
        "d542b1e8dc4bc4dc6f817071b7fd7dc9"),
    "Homo_sapiens.GRCh38.105.gtf.gz": (
        "https://ftp.ensembl.org/pub/release-105/gtf/homo_sapiens/Homo_sapiens.GRCh38.105.gtf.gz", None),
}


def download(url, path):
    if path.exists():
        return
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=60) as response, path.with_suffix(path.suffix + ".partial").open("wb") as out:
                while block := response.read(1024 * 1024):
                    out.write(block)
            path.with_suffix(path.suffix + ".partial").replace(path)
            return
        except Exception:
            if attempt == 2:
                raise
            time.sleep(2)


def write_csv(path, records, fields=None):
    fields = fields or list(dict.fromkeys(key for row in records for key in row))
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(records)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("data_root", type=Path)
    parser.add_argument("--offline", action="store_true")
    args = parser.parse_args()
    out = args.data_root / "molecular_qtl_followup" / "event_annotation"
    raw = out / "raw"
    raw.mkdir(parents=True, exist_ok=True)
    manifests = []
    for name, (url, expected_md5) in SOURCES.items():
        path = raw / name
        if not args.offline:
            download(url, path)
        md5 = hashlib.md5(path.read_bytes()).hexdigest()
        if expected_md5:
            assert md5 == expected_md5, name
        manifests.append({"file": name, "url": url, "bytes": path.stat().st_size,
                          "md5": md5, "expected_md5": expected_md5 or "not supplied",
                          "sha256": hashlib.sha256(path.read_bytes()).hexdigest()})
    write_csv(out / "Annotation_Input_Manifest.csv", manifests)

    events, gff_lines = [], []
    with tarfile.open(raw / "Homo_sapiens.GRCh38.105.txrevise.CAGE25.tar.gz") as archive:
        for member in archive:
            if not member.isfile() or not member.name.endswith(".gff3"):
                continue
            with archive.extractfile(member) as stream:
                for blob in stream:
                    line = blob.decode().rstrip()
                    if GENE not in line:
                        continue
                    gff_lines.append(line)
                    parts = line.split("\t")
                    attrs = dict(field.split("=", 1) for field in parts[8].split(";") if "=" in field)
                    parent = attrs.get("Parent", attrs.get("ID", ""))
                    events.append({"source_member": member.name, "chromosome": parts[0], "feature": parts[2],
                                   "start": int(parts[3]), "end": int(parts[4]), "strand": parts[6],
                                   "event_id": parent, "feature_id": attrs.get("ID", ""),
                                   "attributes": parts[8], "primary_event": parent == EVENT,
                                   "overlaps_rs34536443": int(parts[3]) <= VARIANT_POSITION <= int(parts[4])})
    assert any(row["primary_event"] for row in events)
    (out / "TYK2_frozen_txrevise_events.gff3").write_text("##gff-version 3\n" + "\n".join(gff_lines) + "\n")
    write_csv(out / "TYK2_Txrevise_Event_Features.csv", events)

    with gzip.open(raw / "txrevise_Ensembl_105_phenotype_metadata.tsv.gz", "rt") as stream:
        metadata = [row for row in csv.DictReader(stream, delimiter="\t") if row.get("gene_id") == GENE]
    write_csv(out / "TYK2_Txrevise_Phenotype_Metadata.csv", metadata)

    transcripts, gene_lines = [], []
    with gzip.open(raw / "Homo_sapiens.GRCh38.105.gtf.gz", "rt") as stream:
        for line in stream:
            if f'gene_id "{GENE}"' not in line:
                continue
            gene_lines.append(line.rstrip())
            parts = line.rstrip().split("\t")
            attrs = dict(re.findall(r'(\w+) "([^"]*)"', parts[8]))
            transcripts.append({"chromosome": parts[0], "feature": parts[2], "start": int(parts[3]),
                                "end": int(parts[4]), "strand": parts[6], "phase": parts[7], **attrs})
    (out / "TYK2_Ensembl105.gtf").write_text("\n".join(gene_lines) + "\n")
    write_csv(out / "TYK2_Ensembl105_Transcript_Features.csv", transcripts)
    primary = [row for row in events if row["primary_event"] and row["feature"] == "exon"]
    peers = [row for row in events if row["event_id"].startswith(GENE + ".grp_2.contained.") and row["feature"] == "exon"]
    summaries = []
    for event in sorted({row["event_id"] for row in peers}):
        exons = sorted((row["start"], row["end"]) for row in peers if row["event_id"] == event)
        summaries.append({"event_id": event, "n_exons": len(exons), "exonic_length": sum(b-a+1 for a,b in exons),
                          "intervals_ascending_genomic": ";".join(f"{a}-{b}" for a,b in exons),
                          "primary_event": event == EVENT,
                          "variant_in_event_exon": any(a <= VARIANT_POSITION <= b for a,b in exons)})
    write_csv(out / "TYK2_grp2_Contained_Structure_Summary.csv", summaries)
    audit = {"gene": GENE, "event": EVENT, "event_class": "contained/internal exon structure",
             "annotation_release": "Ensembl105/CAGE25", "strand": primary[0]["strand"],
             "primary_exonic_bases": sum(row["end"]-row["start"]+1 for row in primary),
             "primary_exons": len(primary), "variant_position_GRCh38": VARIANT_POSITION,
             "primary_variant_overlap": any(row["overlaps_rs34536443"] for row in primary),
             "interpretation_limit": "Synthetic event annotation is not a full-length transcript or proof of mediation"}
    (out / "TYK2_Event_Annotation_Audit.json").write_text(json.dumps(audit, indent=2) + "\n")
    print(json.dumps(audit))
    print(out)


if __name__ == "__main__":
    main()
