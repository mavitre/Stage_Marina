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

LiftOn est un outil de lift-over basé sur l'homologie qui intègre des alignements ADN et protéiques pour améliorer la précision de l'annotation à l'échelle du génome, permettant le transfert d'annotations entre espèces relativement distantes. Il s'appuie sur Liftoff et miniprot, et utilise un algorithme de maximisation des protéines pour choisir les meilleurs cadres de lecture ouverts et résoudre les loci géniques chevauchants.


## 3_EGAPx

EGAPx est la version publique du pipeline d'annotation du génome eucaryote du NCBI. Il prend en entrée un fichier FASTA d'assemblage, un taxid de l'organisme, et des données RNA-seq. En fonction du taxid, EGAPx sélectionne automatiquement les ensembles de protéines et les modèles HMM appropriés. Le pipeline utilise miniprot pour aligner les séquences protéiques, STAR pour les lectures RNA-seq courtes, et minimap2 pour les lectures longues, avant de transmettre ces alignements à Gnomon pour la prédiction de gènes

Pour pouvoir lancer EGAPx il est nécessaire d'avoir la bonne version de nextflow. IL n'est pas possible de mettre le nextflow dans pixi il faut donc lancer :
````
NXF_VER=23.10.1 nextflow -version
````
en ligne de commande dans la console du cluster.

Pour pouvoir faire l'annotation il faut remplir un fichier .yaml d'entrée qui est constitué comme ceci:

````
genome: path/to/genome
taxid: trouvable via ce lien : https://www.ncbi.nlm.nih.gov/taxonomy
short_reads: données RNA-seq disponible sur ce lien https://www.ncbi.nlm.nih.gov/sra
 - SRR8506572
 - SRR9005248
cmsearch:
  enabled: false
trnascan:
  enabled: false
````


## 4_helixer

Helixer est un outil ab initio de prédiction de gènes qui produit des modèles géniques précis pour des génomes fongiques, végétaux, vertébrés et invertébrés. Helixer nécessite seulement l'assemblage du génome en fasta, ce qui le rend applicable à une grande diversité d'espèces. Il combine des réseaux de neurones profonds et un modèle de Markov caché pour produire directement des modèles de gènes primaires au format GFF3


## 5_Orthofinder

OrthoFinder est un outil d'inférence d'orthologues conçu pour identifier les gènes orthologues entre plusieurs espèces avec une grande précision.

Script pour recupéré seulement les cds (ou protéines) grâce aux gff produit et leurs fichiers fasta.

Scirpt pour lancer orthofinder

## 6_Post_Orthofinder

### list_orthogroups.py

Les Orthogroups sont trouvé à partir du fichier Orthogroups.GeneCount.tsv dans le dossier Orthogroups.

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


### fetch_orthogroups_sequences.py

Récupère les séquences fasta (dossier "Orthogroup_Sequences" d'OrthoFinder)
correspondant à une liste d'orthogroupes, et les copie dans un dossier de sortie.

Usage :
    python fetch_orthogroups_sequences.py <dossier_resultats_orthofinder> sortie_de_list_orthogroups.txt

    -> Cherche automatiquement le dossier "Orthogroup_Sequences" dans l'arborescence
       (ex: Results_Mmm01/Orthogroup_Sequences/), puis copie les fichiers
       OGxxxxxxx.fa correspondant aux identifiants listés dans orthogroups_1gene.txt
       vers un dossier de sortie.

Options :
    --output D   dossier de sortie où copier les fasta (défaut : orthogroups_sequences_extraites)



## 7_verif_gff

### AGAT

agat_convert_sp_gxf2gxf.pl -g infile.gff [ -o outfile ]

### search_cds_stop_3.py

Pour chaque ARNm (transcrit), fusionne tous ses blocs CDS dans l'ordre
génomique, puis vérifie :
  - que la longueur totale est multiple de 3
  - qu'il n'y a pas plus d'un codon stop dans la séquence traduite
  - que la séquence commence par un codon start (ATG)
  - que la séquence se termine par un codon stop

Se lance avec la commande :
````
python validate_cds.py genome.fasta annotation.gff prefix
````
Pour prefix : donner le nom de l'espèce pour nommer le fichier de sortie

### remove_defective.py

Supprime tous les mRNA défectueux SAUF ceux qui :
  - ont uniquement un problème d'abscence de codon STOP final
  - ET une longueur > 200 nt

````
Usage: python filter_gff_bad_cds.py annotation.gff3 defective_mrna.tsv filtered.gff3
````
defective_mrna.tsv est la sortie de search_cds_stop_3.py
