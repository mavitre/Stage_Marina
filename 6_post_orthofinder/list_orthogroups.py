#!/usr/bin/env python3
"""
Extraction des orthogroupes selon un nombre de gènes par espèce,
DANS CHACUNE des espèces sélectionnées par l'utilisateur.

Usage :
    python list_orthogroups.py Orthogroups.GeneCount.tsv Ae_Bicornis Ae_Comosa Orge

    -> Par défaut (min=max=1) : vérifie que Ae_Bicornis = 1 ET Ae_Comosa = 1 ET Orge = 1
       (orthogroupes "single-copy" stricts).
       Les autres espèces (non listées) ne sont pas prises en compte dans le filtre.

    Pour accepter une plage (ex. entre 1 et 3 gènes par espèce) :
    python list_orthogroups.py Orthogroups.GeneCount.tsv Ae_Bicornis Orge --min 1 --max 3

Options :
    --min N      nombre minimum de gènes par espèce (défaut = 1)
    --max N      nombre maximum de gènes par espèce (défaut = valeur de --min, donc exact par défaut)
    --output F   fichier de sortie (liste des Orthogroup IDs), défaut = orthogroups_1gene.txt
"""

import sys
import argparse
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description="Filtrer les orthogroupes selon un nombre de gènes précis, sur un sous-ensemble d'espèces.")
    parser.add_argument("fichier", help="Fichier de comptage OrthoFinder (Orthogroups.GeneCount.tsv), séparé par tabulations")
    parser.add_argument("especes", nargs="+", help="Noms des espèces à considérer (doivent correspondre exactement aux en-têtes de colonnes)")
    parser.add_argument("--min", type=int, default=1, help="Nombre minimum de gènes par espèce (défaut : 1)")
    parser.add_argument("--max", type=int, default=None, help="Nombre maximum de gènes par espèce ")
    parser.add_argument("--output", default="orthogroups.txt", help="Fichier de sortie pour la liste des Orthogroup IDs")
    args = parser.parse_args()

    # Si --max n'est pas précisé, on cherche une valeur exacte (max = min)
    if args.max is None:
        args.max = args.min

    if args.min > args.max:
        print("Erreur : --min ne peut pas être supérieur à --max")
        sys.exit(1)

    # Lecture du fichier
    df = pd.read_csv(args.fichier, sep="\t")

    # Vérification que les espèces demandées existent bien dans le fichier
    colonnes_dispo = set(df.columns) - {"Orthogroup", "Total"}
    especes_invalides = [e for e in args.especes if e not in colonnes_dispo]
    if especes_invalides:
        print(f"Erreur : ces noms d'espèces n'existent pas dans le fichier : {especes_invalides}")
        print(f"Espèces disponibles : {sorted(colonnes_dispo)}")
        sys.exit(1)

    # Filtrage : orthogroupes où CHAQUE espèce sélectionnée a un nombre de gènes
    # compris entre --min et --max (inclus)
    sous_df = df[args.especes]
    masque = ((sous_df >= args.min) & (sous_df <= args.max)).all(axis=1)
    resultat = df[masque]

    ids = resultat["Orthogroup"].tolist()

    print(f"Espèces sélectionnées : {args.especes}")
    if args.min == args.max:
        critere = f"exactement {args.min} gène(s)"
    else:
        critere = f"entre {args.min} et {args.max} gène(s)"
    print(f"Nombre d'orthogroupes avec {critere} DANS CHACUNE de ces espèces : {len(ids)}")
    print("Identifiants :")
    for i in ids:
        print(f"  {i}")

    # Sauvegarde dans un fichier
    with open(args.output, "w") as f:
        for i in ids:
            f.write(i + "\n")
    print(f"\nListe sauvegardée dans : {args.output}")


if __name__ == "__main__":
    main()