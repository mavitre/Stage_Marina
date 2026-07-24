# Stage Marina

## 1_standardize_gff

Ce code permet de standardizer des gff sur plusieurs critères. 
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

Le seul fichier obligatoire à fournir est le .gff. Il est aussi possible d'ajouter un fichier chromosome.txt en entrée pour ne garder que les chromosomes voulu.
Ce fichier doit se présenter sous la forme suivante:

````
Chr1
Chr2
...
````
Il Faut au préalable prendre en compte la notation des chromosomes dans le gff.

Ce script doit être lancer avec la commande suivante:
````
python standardize_gff.py input.gff [output.gff] [--chromosomes chromosome.txt]
````
Il est conseillé de le lancer en sbatch afin de pouvoir donner une mémoire suffisante.


## 2_LiftOn




## 3_EGAPx

Pour pouvoir lancer EGAPx il est nécessaire d'avoir la bonne version de nextflow. IL n'est pas possible de mettre le nextflow dans pixi il faut donc lancer :
````
NXF_VER=23.10.1 nextflow -version
````
en ligne de commande dans la console du cluster.


## 4_helixer




## 5_Orthofinder



## 6_Post_Orthofinder

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

