# Simplified_scRNAseq
Pour Nicolas Chapelle

***

### En sortie de séquençage, on a des fichiers BCL qu'on veut transformer en FASTQ. Il existe deux méthodes :

## [mkfastq](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/mkfastq)
Nécessite de créer le fichier **index.csv**

## bcl2fastq 
Nécessite de créer le fichier **index_illu.csv**

***

### Puis on utilise les FASTQ pour aligner, démultiplexer les reads, compter l'expression avec cellranger count :

## [cellranger count](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count) (V4,V5) 
requires **library.csv** and **featureindex.csv** when using antibody labelling

***

### Une fois les matrices récupérées, pour Seurat, on utilisera les cellules "pré-filtrées" (alors que pour scDissector, on préferera les cellules "raw")

En sortie de cellranger on a ces outputs à utiliser :

Projet1/    
-⠀⠀|___ outs/  
- - |___ raw_feature_bc_matrix/  
           ⠀⠀⠀⠀⠀⠀⠀⠀⠀|___ matrix.mtx.gz    
⠀⠀⠀⠀⠀⠀⠀⠀⠀           |___ barcodes.tsv.gz   
⠀⠀⠀⠀⠀⠀⠀⠀⠀           |___ features.tsv.gz   
- -  |___ filtered_feature_bc_matrix/  
⠀⠀⠀⠀⠀⠀⠀⠀⠀           |___ matrix.mtx.gz    
⠀⠀⠀⠀⠀⠀⠀⠀⠀           |___ barcodes.tsv.gz   
⠀⠀⠀⠀⠀⠀⠀⠀⠀           |___ features.tsv.gz    

Mon script R "The "PartOne.R" est entièrement automatisé:  
-charge la matrice d'expression  
-seuil pour garder les gènes exprimés dans >3 cells  
-seuil pour garder les cellules avec >200 gènes, <10% mito gènes  
-charger la matrice de comptes HTO et garder les cellules communes à la matrice de compte d'ARNs  
-log normalise par un facteur 10000 et scale les data d'ARN, récupère les top variable genes  
-transforme les data HTO par un log ratio centré  
-demultiplexage des HTO par MULTIseqDemux _( /!\ si les résultats ne sont pas bons, tenter avec HTOdemux)_  
-garde les singlets  
-transforme les data ADT par un log ratio centré  
-enlève les cellules avec un total read count ou unnombre de gènes plus de deux fois la médiane, pour supprimer les doublets potentiels _( /!\ si autre type cellulaire que des PBMC, il est préférable de savoir si il y a une forte hétérogénéité d'expression entre les types cellulaires présents avant d'utiliser ce seuil)_  
-réduction de dimensions en PCA  
-decision sur le nombre de PC à utiliser pour la suite, projection en UMAP  
-clustree dessine un arbre du nombre de cluster par resolutions, pour aider à décider du nombre de cluster optimal  

