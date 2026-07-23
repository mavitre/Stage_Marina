#!/usr/bin/env python3
"""
standardize_gff.py — Standardisation et filtrage de fichiers GFF3

Règles appliquées, dans cet ordre :

  0. Chromosomes (optionnel, si --chromosomes est fourni)
     Ne garde que les features dont la colonne 1 (seqid) figure dans le
     fichier .txt fourni (une valeur par ligne). 

  1. Type de feature (colonne 3)
     Ne conserve que : gene, mRNA, CDS, exon, five_prime_UTR, three_prime_UTR

  2. Confidence (High / Low)
     Détection automatique de l'attribut utilisé : "confidence" ou "primconf"
     - mRNA Low                          → supprimé
     - gene Low SANS aucun mRNA High     → supprimé (gène + tous ses descendants)
     - gene Low AVEC au moins un mRNA High → conservé (seuls ses mRNA Low sont retirés)
     - gene High qui se retrouve sans aucun mRNA restant après filtrage → supprimé

  3. Ensembl_canonical
     Si le tag "tag=Ensembl_canonical" est présent quelque part dans le fichier,
     ne garde que les mRNA canoniques.

  4. Note de modification de séquence
     Supprime toute feature (et ses descendants) portant la note :
     "Note=The sequence of the model RefSeq protein was modified"

  5. Nettoyage final des gènes orphelins
     Après tous les filtres ci-dessus, supprime tout gene qui n'a plus
     aucun mRNA enfant restant.

  6. Propagation récursive
     Toute suppression (gene ou mRNA) entraîne la suppression de TOUS ses descendants.

Usage :
    python standardize_gff.py input.gff [output.gff] [--chromosomes chrom_list.txt]

    Si output.gff est omis : input.standardized.gff
    --chromosomes / -c : fichier .txt listant (une par ligne) les valeurs de
                          la colonne seqid (chromosome/scaffold) à conserver.
                          Si absent, aucun filtre sur les chromosomes n'est appliqué.
"""

import sys
import os
import argparse


KEPT_FEATURES = {"gene", "mRNA", "CDS", "exon", "five_prime_UTR", "three_prime_UTR"}

# Noms d'attributs pouvant porter l'information de confiance, par ordre de priorité
CONFIDENCE_ATTR_CANDIDATES = ["confidence", "primconf"]

LOW_VALUES = {"low", "low_confidence", "lc"}
HIGH_VALUES = {"high", "high_confidence", "hc"}

# Note signalant une modification de séquence du modèle RefSeq → feature à exclure
MODIFIED_SEQUENCE_NOTE = "the sequence of the model refseq protein was modified"

# ──────────────────────────────────────────────────────────────────────────
# Parsing
# ──────────────────────────────────────────────────────────────────────────

def parse_attributes(attr_string):
    attrs = {}
    for part in attr_string.strip().split(";"):
        part = part.strip()
        if not part:
            continue
        if "=" in part:
            key, _, value = part.partition("=")
            attrs[key.strip()] = value.strip()
        else:
            attrs[part] = ""
    return attrs

def get_parents(attrs):
    """Retourne la liste des Parent IDs (gère les Parent multiples séparés par virgule)."""
    parent = attrs.get("Parent", "")
    return [p.strip() for p in parent.split(",") if p.strip()]

def has_canonical_tag(attrs):
    tags = attrs.get("tag", "")
    return "Ensembl_canonical" in [t.strip() for t in tags.split(",")]

def has_modified_sequence_note(attrs):
    """Détecte la note 'The sequence of the model RefSeq protein was modified'
    (insensible à la casse, et tolère un contenu plus large dans Note=)."""
    note = attrs.get("Note", "")
    return MODIFIED_SEQUENCE_NOTE in note.lower()

def get_confidence(attrs):
    """Retourne 'high', 'low' ou None selon les attributs disponibles."""
    for attr_name in CONFIDENCE_ATTR_CANDIDATES:
        if attr_name in attrs:
            val = attrs[attr_name].strip().lower()
            if val in HIGH_VALUES:
                return "high"
            if val in LOW_VALUES:
                return "low"
    return None

# ──────────────────────────────────────────────────────────────────────────
# Chargement
# ──────────────────────────────────────────────────────────────────────────

def load_gff(filepath):
    """
    Retourne :
      headers : lignes de commentaires/directives/malformées, à réécrire telles quelles
      records : liste de dicts {line, fields, attrs, id, parents, feature, uid}

    """
    headers = []
    records = []
    uid_counter = 0
    with open(filepath, "r", encoding="utf-8") as fh:
        for line in fh:
            raw = line.rstrip("\n")
            if raw.startswith("#") or raw.strip() == "":
                headers.append(raw)
                continue
            fields = raw.split("\t")
            if len(fields) < 9:
                headers.append(raw)
                continue
            attrs = parse_attributes(fields[8])
            uid_counter += 1
            records.append({
                "line": raw,
                "fields": fields,
                "attrs": attrs,
                "id": attrs.get("ID", ""),
                "parents": get_parents(attrs),
                "feature": fields[2],
                "seqid": fields[0],
                "uid": uid_counter,
            })
    return headers, records

def load_chromosome_list(filepath):
    """Charge un fichier .txt (une valeur par ligne) listant les chromosomes à garder."""
    keep = set()
    with open(filepath, "r", encoding="utf-8") as fh:
        for line in fh:
            val = line.strip()
            if val and not val.startswith("#"):
                keep.add(val)
    return keep

# ──────────────────────────────────────────────────────────────────────────
# Étape 0 : filtre par chromosome (optionnel)
# ──────────────────────────────────────────────────────────────────────────

def filter_by_chromosome(records, chrom_list):
    """Ne garde que les features dont seqid (colonne 1) figure dans chrom_list."""
    if not chrom_list:
        return records, {"applied": False}
    filtered = [r for r in records if r["seqid"] in chrom_list]
    return filtered, {"applied": True, "removed": len(records) - len(filtered)}

# ──────────────────────────────────────────────────────────────────────────
# Étape 1 : filtre par type de feature
# ──────────────────────────────────────────────────────────────────────────

def filter_by_feature_type(records):
    return [r for r in records if r["feature"] in KEPT_FEATURES]

# ──────────────────────────────────────────────────────────────────────────
# Étape 2 : filtre confidence (High / Low)
# ──────────────────────────────────────────────────────────────────────────

def filter_by_confidence(records):
    """
    Retourne (records_filtrés, stats_dict).
    """
    detected = any(get_confidence(r["attrs"]) is not None for r in records)
    if not detected:
        return records, {"applied": False}

    # gene_id / mrna_id (attribut GFF ID=) → uid de la ligne correspondante
    # (gene et mRNA ont normalement un ID unique et non vide)
    gene_uid_by_id = {r["id"]: r["uid"] for r in records if r["feature"] == "gene" and r["id"]}
    mrna_uid_by_id = {r["id"]: r["uid"] for r in records if r["feature"] == "mRNA" and r["id"]}

    gene_conf = {r["id"]: get_confidence(r["attrs"]) for r in records if r["feature"] == "gene" and r["id"]}
    mrna_conf = {r["id"]: get_confidence(r["attrs"]) for r in records if r["feature"] == "mRNA" and r["id"]}

    # gene_id → liste des mRNA_id enfants (résolution par ID GFF, comme le veut le format)
    mrna_children_of_gene = {}
    for r in records:
        if r["feature"] == "mRNA":
            for p in r["parents"]:
                mrna_children_of_gene.setdefault(p, []).append(r["id"])

    to_remove_ids = set()  # ensemble d'ID GFF (gene_id / mrna_id) marqués supprimés

    # gènes Low sans aucun mRNA High → supprimés
    for gene_id, conf in gene_conf.items():
        if conf == "low":
            mrna_ch = mrna_children_of_gene.get(gene_id, [])
            has_high_mrna = any(mrna_conf.get(m) == "high" for m in mrna_ch)
            if not has_high_mrna:
                to_remove_ids.add(gene_id)

    # tous les mRNA Low → supprimés
    for mrna_id, conf in mrna_conf.items():
        if conf == "low":
            to_remove_ids.add(mrna_id)

    # gènes qui se retrouvent sans aucun mRNA restant → supprimés
    for gene_id in gene_conf:
        if gene_id in to_remove_ids:
            continue
        mrna_ch = mrna_children_of_gene.get(gene_id, [])
        if mrna_ch and all(m in to_remove_ids for m in mrna_ch):
            to_remove_ids.add(gene_id)

    filtered, removed_uids = remove_by_ancestor_ids(records, to_remove_ids)

    stats = {
        "applied": True,
        "low_genes_removed": sum(1 for g in to_remove_ids if gene_conf.get(g) == "low"),
        "high_genes_removed": sum(1 for g in to_remove_ids if gene_conf.get(g) == "high"),
        "low_mrna_removed": sum(1 for m in to_remove_ids if mrna_conf.get(m) == "low"),
    }
    return filtered, stats

# ──────────────────────────────────────────────────────────────────────────
# Étape 3 : filtre Ensembl_canonical
# ──────────────────────────────────────────────────────────────────────────

def filter_by_canonical(records):
    """
    Si tag=Ensembl_canonical détecté sur au moins un mRNA, ne garde que :
      - les mRNA canoniques
      - les gènes parents ayant encore au moins un mRNA canonique
      - tous les descendants des mRNA canoniques
    """
    canonical_mrna_ids = {
        r["id"] for r in records
        if r["feature"] == "mRNA" and has_canonical_tag(r["attrs"]) and r["id"]
    }
    if not canonical_mrna_ids:
        return records, {"applied": False}

    genes_to_keep = set()
    for r in records:
        if r["feature"] == "mRNA" and r["id"] in canonical_mrna_ids:
            genes_to_keep.update(r["parents"])

    to_remove_ids = {
        r["id"] for r in records
        if r["feature"] == "mRNA" and r["id"] not in canonical_mrna_ids and r["id"]
    }
    for r in records:
        if r["feature"] == "gene" and r["id"] not in genes_to_keep and r["id"]:
            to_remove_ids.add(r["id"])

    filtered, _ = remove_by_ancestor_ids(records, to_remove_ids)

    stats = {"applied": True, "canonical_mrna_kept": len(canonical_mrna_ids)}
    return filtered, stats

# ──────────────────────────────────────────────────────────────────────────
# Étape 4 : filtre Note "sequence of the model RefSeq protein was modified"
# ──────────────────────────────────────────────────────────────────────────

def filter_by_modified_note(records):
    """
    Supprime toute feature portant cette Note, ainsi que tous ses descendants.
    """
    direct_uids_to_remove = {
        r["uid"] for r in records if has_modified_sequence_note(r["attrs"])
    }
    ids_to_remove = {
        r["id"] for r in records
        if has_modified_sequence_note(r["attrs"]) and r["id"]
    }
    if not direct_uids_to_remove:
        return records, {"applied": False, "removed": 0}

    filtered, removed_uids = remove_by_ancestor_ids(
        records, ids_to_remove, extra_uids=direct_uids_to_remove
    )

    stats = {"applied": True, "removed": len(removed_uids)}
    return filtered, stats


# ──────────────────────────────────────────────────────────────────────────
# Étape 5 : nettoyage final des gènes orphelins (sans aucun mRNA)
# ──────────────────────────────────────────────────────────────────────────

def remove_orphan_genes(records):
    """
    Supprime tout gene qui, après tous les filtres précédents, n'a plus
    aucun mRNA enfant.
    """
    mrna_ids_by_parent = {}
    for r in records:
        if r["feature"] == "mRNA":
            for p in r["parents"]:
                mrna_ids_by_parent.setdefault(p, []).append(r["id"])

    orphan_gene_ids = {
        r["id"] for r in records
        if r["feature"] == "gene" and r["id"]
        and not mrna_ids_by_parent.get(r["id"])
    }

    if not orphan_gene_ids:
        return records, {"applied": False, "orphan_genes_removed": 0}

    filtered, removed_uids = remove_by_ancestor_ids(records, orphan_gene_ids)
    removed_total = len(orphan_gene_ids)

    return filtered, {"applied": removed_total > 0, "orphan_genes_removed": removed_total}


# ──────────────────────────────────────────────────────────────────────────
# Suppression robuste par ID GFF ancêtre (gère ID vides/partagés correctement)
# ──────────────────────────────────────────────────────────────────────────

def remove_by_ancestor_ids(records, ancestor_ids, extra_uids=None):
    """
    Supprime toutes les lignes dont l'ID GFF (ou un de ses ancêtres via Parent=)
    figure dans `ancestor_ids`, en propageant récursivement à toute profondeur.

    Cette fonction travaille en UID interne (un par ligne), jamais en ID GFF brut,
    pour éviter toute collision quand des lignes (typiquement exon/UTR) ont un
    ID= vide ou partagé (typiquement CDS multi-exons partageant le même ID).

    Retourne (records_filtrés, set_des_uid_supprimés).
    """
    ancestor_ids = {a for a in ancestor_ids if a}  # ignore les ID vides en entrée

    # id GFF parent → liste d'uid des lignes enfants directes
    children_uids_by_parent_id = {}
    for r in records:
        for p in r["parents"]:
            children_uids_by_parent_id.setdefault(p, []).append(r["uid"])

    # uid → id GFF de la ligne (pour pouvoir continuer la propagation depuis un enfant)
    id_by_uid = {r["uid"]: r["id"] for r in records}

    to_remove_uids = set(extra_uids) if extra_uids else set()

    # initialiser avec toutes les lignes dont l'ID figure dans ancestor_ids
    for r in records:
        if r["id"] and r["id"] in ancestor_ids:
            to_remove_uids.add(r["uid"])

    # propagation récursive via les ID GFF (un enfant peut lui-même être parent)
    stack = [id_by_uid[u] for u in to_remove_uids if id_by_uid.get(u)]
    seen_ids = set(stack)
    while stack:
        current_id = stack.pop()
        for child_uid in children_uids_by_parent_id.get(current_id, []):
            if child_uid not in to_remove_uids:
                to_remove_uids.add(child_uid)
                child_id = id_by_uid.get(child_uid)
                if child_id and child_id not in seen_ids:
                    seen_ids.add(child_id)
                    stack.append(child_id)

    filtered = [r for r in records if r["uid"] not in to_remove_uids]
    return filtered, to_remove_uids

# ──────────────────────────────────────────────────────────────────────────
# Écriture
# ──────────────────────────────────────────────────────────────────────────

def write_output(output_path, headers, records):
    with open(output_path, "w", encoding="utf-8") as fh:
        for h in headers:
            fh.write(h + "\n")
        for r in records:
            fh.write(r["line"] + "\n")

# ──────────────────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────────────────

def count_genes(records):
    return sum(1 for r in records if r["feature"] == "gene")


def main():
    parser = argparse.ArgumentParser(
        description="Standardisation et filtrage de fichiers GFF3",
        add_help=True,
    )
    parser.add_argument("input", help="Fichier GFF3 d'entrée")
    parser.add_argument("output", nargs="?", default=None,
                         help="Fichier GFF3 de sortie (optionnel)")
    parser.add_argument("-c", "--chromosomes", default=None,
                         help="Fichier .txt listant les chromosomes/seqid à conserver "
                              "(un par ligne). Si omis, aucun filtre chromosome n'est appliqué.")
    args = parser.parse_args()

    input_path = args.input
    if args.output:
        output_path = args.output
    else:
        base, ext = os.path.splitext(input_path)
        output_path = base + ".standardized" + (ext if ext else ".gff")

    if not os.path.isfile(input_path):
        print(f"[ERREUR] Fichier introuvable : {input_path}", file=sys.stderr)
        sys.exit(1)

    chrom_list = None
    if args.chromosomes:
        if not os.path.isfile(args.chromosomes):
            print(f"[ERREUR] Fichier chromosomes introuvable : {args.chromosomes}", file=sys.stderr)
            sys.exit(1)
        chrom_list = load_chromosome_list(args.chromosomes)
        print(f"Liste de chromosomes chargée : {len(chrom_list)} chromosome.s")

    print(f"Chargement de : {input_path}")
    headers, records = load_gff(input_path)
    print(f"Gènes initiaux : {count_genes(records)}")

    # Étape 0 — chromosomes
    records, chrom_stats = filter_by_chromosome(records, chrom_list)
    if chrom_stats["applied"]:
        print(f"Filtre chromosomes appliqué")
        print(f"Gènes après filtre chromosomes : {count_genes(records)}")
    else:
        print(f"Pas de filtre chromosomes (--chromosomes non fourni)")
    
    # Étape 1 — type de feature
    records = filter_by_feature_type(records)

    # Étape 2 — confidence High/Low
    records, conf_stats = filter_by_confidence(records)
    if conf_stats["applied"]:
        print(f"Filtre confidence appliqué (High/Low détecté)")
        print(f"       Gènes Low supprimés          : {conf_stats['low_genes_removed']}")
        print(f"       Gènes High supprimés (vides) : {conf_stats['high_genes_removed']}")
        print(f"       mRNA Low supprimés           : {conf_stats['low_mrna_removed']}")
        print(f"Gènes après filtre confidence : {count_genes(records)}")
    else:
        print(f"Pas de confidence/primconf détecté → filtre ignoré")
    
    # Étape 3 — Ensembl_canonical
    records, canon_stats = filter_by_canonical(records)
    if canon_stats["applied"]:
        print(f"Filtre Ensembl_canonical appliqué")
        print(f"       mRNA canoniques conservés : {canon_stats['canonical_mrna_kept']}")
    else:
        print(f"Pas de tag Ensembl_canonical détecté")
    

    # Étape 4 — Note "sequence of the model RefSeq protein was modified"
    records, note_stats = filter_by_modified_note(records)
    if note_stats["applied"]:
        print(f"Filtre Note 'modified RefSeq protein' appliqué")
        print(f"Gènes après filtre Note : {count_genes(records)}")
    else:
        print(f"Pas de Note 'modified RefSeq protein' détectée")
    
    # Étape 5 — nettoyage des gènes orphelins (sans mRNA)
    records, orphan_stats = remove_orphan_genes(records)
    if orphan_stats["applied"]:
        print(f"Nettoyage gènes sans mRNA appliqué")
        print(f"       Gènes sans mRNA supprimés : {orphan_stats['orphan_genes_removed']}")
    else:
        print(f"Pas de gènes sans descendants")

    write_output(output_path, headers, records)
    print(f"Gènes finaux : {count_genes(records)}")
    print(f"Fichier écrit : {output_path}")


if __name__ == "__main__":
    main()