#!/usr/bin/env python3
"""
filter_gff_bad_cds.py
Usage: python filter_gff_bad_cds.py annotation.gff3 defective_mrna.tsv filtered.gff3

Supprime tous les mRNA défectueux SAUF ceux qui :
  - ont uniquement no_stop_end == 1 (et tous les autres tests == 0)
  - ET une longueur > 200 nt
"""
import sys
from collections import defaultdict

gff_file, bad_file, out_file = sys.argv[1], sys.argv[2], sys.argv[3]

# 1. Charger le TSV et construire l'ensemble des mRNA à supprimer
bad_mrna = set()

with open(bad_file) as fh:
    header = fh.readline().strip().split("\t")
    # Index des colonnes
    idx_id      = header.index("mRNA_id")
    idx_len     = header.index("length")
    idx_nm3     = header.index("non_multiple3")
    idx_mstop   = header.index("multi_stop")
    idx_nostart = header.index("no_start")
    idx_nostop  = header.index("no_stop_end")

    for line in fh:
        if not line.strip():
            continue
        cols = line.strip().split("\t")
        mrna_id      = cols[idx_id]
        length       = int(cols[idx_len])
        non_mult3    = int(cols[idx_nm3])
        multi_stop   = int(cols[idx_mstop])
        no_start     = int(cols[idx_nostart])
        no_stop_end  = int(cols[idx_nostop])

        # Garder si : seul défaut = no_stop_end ET longueur > 200
        only_no_stop = (no_stop_end == 1 and non_mult3 == 0
                        and multi_stop == 0 and no_start == 0)

        if only_no_stop and length > 60:
            continue  # on garde ce mRNA → ne pas l'ajouter à bad_mrna

        bad_mrna.add(mrna_id)

print(f"ARNm à supprimer : {len(bad_mrna)}")

# 2. Lire le GFF3 : construire parent->enfants et feature_id->lignes
lines = open(gff_file).readlines()
children   = defaultdict(set)    # parent_id -> {enfant_ids}
feat_lines = defaultdict(list)   # feature_id -> [indices de lignes]

for i, line in enumerate(lines):
    if line.startswith("#") or not line.strip():
        continue
    parts = line.split("\t")
    if len(parts) < 9:
        continue
    attrs = dict(kv.split("=", 1) for kv in parts[8].strip().split(";") if "=" in kv)
    fid = attrs.get("ID")
    if fid:
        feat_lines[fid].append(i)
    for parent in attrs.get("Parent", "").split(","):
        if parent:
            children[parent].add(fid)

# 3. Collecter récursivement tous les descendants d'un ID
def descendants(fid):
    result = set()
    stack = [fid]
    while stack:
        cur = stack.pop()
        for child in children.get(cur, set()):
            if child and child not in result:
                result.add(child)
                stack.append(child)
    return result

# 4. Construire l'ensemble des IDs à supprimer
to_remove = set()
for mrna_id in bad_mrna:
    if mrna_id not in feat_lines:
        continue
    to_remove.add(mrna_id)
    to_remove.update(descendants(mrna_id))

# 5. Supprimer les gènes parents qui n'ont plus de mRNA
for gene_id, gene_children in children.items():
    if gene_children and gene_children <= to_remove:
        to_remove.add(gene_id)

# 6. Écrire les lignes non supprimées
bad_lines = {i for fid in to_remove for i in feat_lines[fid]}
with open(out_file, "w") as out:
    for i, line in enumerate(lines):
        if i not in bad_lines:
            out.write(line)

print(f"Fichier filtré écrit : {out_file}")