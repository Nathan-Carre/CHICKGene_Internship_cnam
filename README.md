# CHIKGene — Differential Expression Analysis (DEG)

## Introduction
Ce dépôt décrit en détail la démarche à suivre pour reproduire les analyses de gènes différentiellement exprimés (DEG) menées dans le cadre du stage de Master 2 sur la cohorte CHIKGene. Tous les scripts et étapes nécessaires y sont présentés afin de permettre une réexécution transparente de l’analyse.
Pour des raisons de confidentialité et de protection des données patients, les fichiers de comptage et les fichiers phénotypes (données cliniques) ne peuvent pas être publiés ici. Le dépôt ne contient donc que les scripts et instructions permettant de les utiliser une fois les données disponibles en interne.


## Contexte
Ce projet regroupe les analyses de gènes différentiellement exprimés (DEG) réalisées dans le cadre du stage de Master 2 sur la cohorte CHIKGene (données PBMCs RNA-seq). Les analyses visent à identifier des signatures transcriptomiques associées à différents symptômes chroniques post-infection Chikungunya (DMS, FCI, cardiopathies, diabète, double-atteinte...)


## Données d’entrée
- **Matrice de comptage** (ex : `CountsChikgeneFinal133.txt` ou 'filtered_cut_219_CountsChikgene.txt.zip' qui doit être dézippé au préalable)
- **Phénotypes** : fichiers CSV par analyse avec au minimum `sex`, `tecnica`, la variable de **condition** ciblée, et éventuellement d’autres covariables (CHIKR, tabac, statut pré/post infection par le Chikungunya).


## Scripts d’analyse (6 comparaisons principales)

Chaque script suit le même schéma : import des données, design DESeq2 avec ajustement (sexe, technique, + ou - la(s) covariable(s) spécifique(s) ), exécution du modèle, export des résultats et matrices normalisées/log-normalisées. Les **Volcano plots** sont optionnels.

1. **01_DMS_CC_noCHIKR_yes_vs_no.R**  
   - Comparaison : CC DMS oui vs non (exclusion CHIKR)
   - Sorties : `resBatchEffectSex_Tecnica_type_symptôme.csv`, `norm_counts_*.csv`, `log_norm_counts_*.csv`...

2. **02_FCI_CC_noCHIKR_yes_vs_no.R**  
   - Comparaison : CC FCI oui vs non (exclusion CHIKR)
   - Sorties : `resBatchEffectSex_Tecnica_type_symptôme.csv`, `norm_counts_*.csv`, `log_norm_counts_*.csv`...
     
3. **04_Cardio_CC_post_SMOKE_yes_vs_no.R**  
   - Comparaison : fumeurs vs non-fumeurs (chez CC avec cardiopathie post-infection)
   - Sorties : `resBatchEffectSex_Tecnica_type_symptôme.csv`, `norm_counts_*.csv`, `log_norm_counts_*.csv`...

4. **05_Diabete_CC_yes_vs_no.R**  
   - Comparaison : diabète oui vs non (chez CC)
   - Sorties : `resBatchEffectSex_Tecnica_type_symptôme.csv`, `norm_counts_*.csv`, `log_norm_counts_*.csv`...

5. **06_Diabete_CC_post_vs_pre.R**  
   - Comparaison : diabète post-infection vs pré-infection (chez CC)
   - Sorties : `resBatchEffectSex_Tecnica_type_symptôme.csv`, `norm_counts_*.csv`, `log_norm_counts_*.csv`...
  
6. **06_(Diabete_Cardio)_yes_vs_noDiabete_noCardio.R**  
   - Comparaison : diabète_et_cardio (les deux à la fois) vs nonDiabète_et_nonCardio (ni l'un ni l'autre) (chez CC)
   - Sorties : `resBatchEffectSex_Tecnica_type_symptôme.csv`, `norm_counts_*.csv`, `log_norm_counts_*.csv`...


## Tri des résultats
- **Sort_Result_DEG.R** : trie chaque fichier `resBatchEffectSex_*_.csv` par p-value croissante et génère un fichier `SORTED_pvalue_*_.csv`.
- À exécuter systématiquement après chaque analyse DEG, si l'on veut ensuite utilisé les gene_name par ordre croissant de p-value dans l'analyse Metascape


## Analyse fonctionnelle et annotation
1. **Metascape** : enrichissement fonctionnel (GO, KEGG, Reactome, WikiPathways, DisGeNET).
2. **DisGeNET** : associations gènes–maladies.
3. **GeneCards** : exploration des gènes clés (expression, fonctions, maladies associées).
4. **NCBI (Gene, PubMed)** : validation par articles.

Les listes de gènes up/down (pvalue < 0.05, |log2FC| > 0.7) issues des fichiers triés servent d’entrée pour ces outils. On sélectionne tous les gene_name (première colonne) par ordre croissant de p-value strictement inférieures à 0.05. Puis on copie/colle cette liste dans l'entrée Metascape, avec comme sélection de 'Input as species:' l'espèce humaine soit 'H. Sapiens (6)'


## Workflow type
1. Lancer le script d’analyse (ex : `CC_DMS_25_06_ChikR.R`) avec les packages installés au préalable
2. (Optionnel) Générer volcano plots et/ou MA-plots
3. Lancer le tri : `Sort_Result_DEG.R` sur la sortie correspondante
5. Exporter les listes gènes up/down selon le critère de p-value pour les mettre en entrées dans Metascape
6. Faire l’enrichissement fonctionnel (Metascape, DisGeNET) et annoter gènes d’intérêt via GeneCards et NCBI


## Conseils
- Bien installer tous les packages nécessaire et vérifier leur bon fonctionnement
- Vérifier le sens du contraste (case vs control) avant interprétation.
- Respecter la construction des matrices (0/1/2) pour exclure les groupes non pertinents dans les design.
