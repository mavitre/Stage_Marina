#!/usr/bin/env python3
"""
Récupère les séquences fasta (dossier "Orthogroup_Sequences" d'OrthoFinder)
correspondant à une liste d'orthogroupes, et les copie dans un dossier de sortie.

Usage :
    python fetch_orthogroup_sequences.py <dossier_resultats_orthofinder> orthogroups_1gene.txt

    -> Cherche automatiquement le dossier "Orthogroup_Sequences" dans l'arborescence
       (ex: Results_Mmm01/Orthogroup_Sequences/), puis copie les fichiers
       OGxxxxxxx.fa correspondant aux identifiants listés dans orthogroups_1gene.txt
       vers un dossier de sortie.

Options :
    --output D   dossier de sortie où copier les fasta (défaut : orthogroups_sequences_extraites)
"""

import sys
import argparse
import shutil
from pathlib import Path


def trouver_dossier_sequences(racine: Path) -> Path:
    """Cherche un dossier nommé 'Orthogroup_Sequences' n'importe où sous 'racine'."""
    candidats = [p for p in racine.rglob("*") if p.is_dir() and p.name == "Orthogroup_Sequences"]
    if not candidats:
        return None
    # S'il y en a plusieurs (plusieurs runs), on prend le plus récent
    return max(candidats, key=lambda p: p.stat().st_mtime)


def lire_liste_ids(fichier_liste: Path) -> list:
    with open(fichier_liste) as f:
        return [ligne.strip() for ligne in f if ligne.strip()]


def main():
    parser = argparse.ArgumentParser(description="Extraire les séquences fasta d'une liste d'orthogroupes depuis l'arborescence OrthoFinder.")
    parser.add_argument("dossier_orthofinder", help="Dossier racine des résultats OrthoFinder (ex: Results_Mmm01, ou son dossier parent)")
    parser.add_argument("liste_ids", help="Fichier texte contenant les identifiants d'orthogroupes (un par ligne), ex: orthogroups_1gene.txt")
    parser.add_argument("--output", default="orthogroups_sequences_extraites", help="Dossier de sortie pour les fasta copiés")
    args = parser.parse_args()

    racine = Path(args.dossier_orthofinder)
    if not racine.exists():
        print(f"Erreur : le dossier '{racine}' n'existe pas.")
        sys.exit(1)

    dossier_seq = trouver_dossier_sequences(racine)
    if dossier_seq is None:
        print(f"Erreur : aucun dossier 'Orthogroup_Sequences' trouvé sous '{racine}'.")
        print("Vérifiez que vous pointez bien vers le dossier de résultats OrthoFinder (ou un de ses parents).")
        sys.exit(1)

    print(f"Dossier de séquences trouvé : {dossier_seq}")

    ids = lire_liste_ids(Path(args.liste_ids))
    print(f"{len(ids)} orthogroupe(s) à récupérer depuis : {args.liste_ids}")

    dossier_sortie = Path(args.output)
    dossier_sortie.mkdir(parents=True, exist_ok=True)

    trouves = []
    manquants = []

    for og_id in ids:
        source = dossier_seq / f"{og_id}.fa"
        if not source.exists():
            # essai avec extension .fasta au cas où
            source_alt = dossier_seq / f"{og_id}.fasta"
            if source_alt.exists():
                source = source_alt
            else:
                manquants.append(og_id)
                continue
        shutil.copy2(source, dossier_sortie / source.name)
        trouves.append(og_id)

    print(f"\n{len(trouves)} fichier(s) copié(s) vers : {dossier_sortie}")
    if manquants:
        print(f"{len(manquants)} orthogroupe(s) introuvable(s) dans '{dossier_seq}' :")
        for m in manquants:
            print(f"  {m}")


if __name__ == "__main__":
    main()