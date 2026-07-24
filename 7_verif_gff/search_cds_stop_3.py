#!/usr/bin/env python3
"""
validate_cds.py
 
Usage: python validate_cds.py genome.fasta annotation.gff prefix
 
Pour chaque ARNm (transcrit), fusionne tous ses blocs CDS dans l'ordre
génomique, puis vérifie :
  - que la longueur totale est multiple de 3
  - qu'il n'y a pas plus d'un codon stop dans la séquence traduite
  - que la séquence commence par un codon start (ATG)
  - que la séquence se termine par un codon stop
 
Dépendance: pip install biopython
"""
import sys
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq

STOP_CODONS = {"TAA", "TAG", "TGA"}
START_CODON = "ATG"


def load_gff(gff_path):
    mrna = defaultdict(lambda: {"seqid": None, "strand": None, "blocks": []})

    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2].upper() != "CDS":
                continue

            attrs = dict(
                f.split("=", 1)
                for f in cols[8].strip().split(";")
                if "=" in f
            )

            parent = attrs.get("Parent", attrs.get("ID", f"{cols[0]}_{cols[3]}_{cols[4]}"))
            parent = parent.split(",")[0]

            seqid  = cols[0]
            strand = cols[6]
            start  = int(cols[3]) - 1
            end    = int(cols[4])

            mrna[parent]["seqid"]  = seqid
            mrna[parent]["strand"] = strand
            mrna[parent]["blocks"].append((start, end))

    for rec in mrna.values():
        rec["blocks"].sort()

    return mrna


def get_merged_cds_seq(rec, genome):
    chrom_seq = genome[rec["seqid"]].seq
    merged = Seq("".join(str(chrom_seq[s:e]) for s, e in rec["blocks"]))

    if rec["strand"] == "-":
        merged = merged.reverse_complement()

    return merged


def count_stops(seq):
    return sum(
        str(seq[i:i+3]).upper() in STOP_CODONS
        for i in range(0, len(seq) - 2, 3)
    )


def main():
    if len(sys.argv) != 4:
        sys.exit("Usage: python validate_cds.py genome.fasta annotation.gff prefix")

    fasta_path, gff_path, prefix = sys.argv[1], sys.argv[2], sys.argv[3]

    print("Chargement du génome…")
    genome = SeqIO.index(fasta_path, "fasta")

    print("Lecture du GFF…")
    mrna_map = load_gff(gff_path)
    print(f"  → {len(mrna_map)} ARNm avec au moins un bloc CDS")

    defective = {}

    for mrna_id, rec in mrna_map.items():
        if rec["seqid"] not in genome:
            print(f"  [WARN] séquence absente du FASTA : {rec['seqid']} (ARNm {mrna_id})")
            continue

        seq       = get_merged_cds_seq(rec, genome)
        total_len = len(seq)

        non_multiple3 = int(total_len % 3 != 0)
        multi_stop    = int(count_stops(seq) > 1)
        no_start      = int(str(seq[:3]).upper() != START_CODON)
        no_stop_end   = int(str(seq[-3:]).upper() not in STOP_CODONS)

        if any([non_multiple3, multi_stop, no_start, no_stop_end]):
            defective[mrna_id] = {
                "seqid":         rec["seqid"],
                "strand":        rec["strand"],
                "length":        total_len,
                "non_multiple3": non_multiple3,
                "multi_stop":    multi_stop,
                "no_start":      no_start,
                "no_stop_end":   no_stop_end,
            }

    # --- Écriture du fichier unique ---
    out_path = f"{prefix}_defective_mrna.tsv"
    with open(out_path, "w") as fh:
        fh.write("mRNA_id\tseqid\tstrand\tlength\tnon_multiple3\tmulti_stop\tno_start\tno_stop_end\n")
        for mrna_id, info in defective.items():
            fh.write(
                f"{mrna_id}\t"
                f"{info['seqid']}\t"
                f"{info['strand']}\t"
                f"{info['length']}\t"
                f"{info['non_multiple3']}\t"
                f"{info['multi_stop']}\t"
                f"{info['no_start']}\t"
                f"{info['no_stop_end']}\n"
            )

    print(f"\nTotal ARNm défectueux : {len(defective)}")
    print(f"Fichier de sortie     : {out_path}")

    print("\nRésumé par type de défaut :")
    print(f"  non_multiple3 : {sum(v['non_multiple3'] for v in defective.values())}")
    print(f"  multi_stop    : {sum(v['multi_stop']    for v in defective.values())}")
    print(f"  no_start      : {sum(v['no_start']      for v in defective.values())}")
    print(f"  no_stop_end   : {sum(v['no_stop_end']   for v in defective.values())}")


if __name__ == "__main__":
    main()
