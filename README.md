# Human glial progenitors transplanted into Huntington Disease mice substantially rescue neuronal gene expression, dendritic structure and behavior 

## **Publication**

Villanueva et al., 2025

## **Abstract**
Neonatal glial replacement can delay disease progression in mouse models of Huntington’s disease (HD), in which glial dysfunction is a prominent feature. Here we asked if transplanting healthy human glial progenitors (hGPCs) into adult R6/2 HD mice might ameliorate phenotype, while concurrently investigating the diseased host’s neuronal response to wild-type glial engraftment. We found that the introduction of hGPCs into the striata of adult R6/2 HD mice indeed delayed their motor and cognitive decline, while extending survival. Single nucleus RNA sequencing revealed that whereas neuronal genes associated with synaptogenesis and synaptic structure were differentially down-regulated in R6/2 striatal neurons, the transcription of these genes was relatively rescued by wild-type glial engraftment. Retrograde rabies labeling of striatal MSNs then revealed that while both dendritic complexity and dendritic spine density were deficient in R6/2 mice, each was largely restored in hGPC-engrafted R6/2 mice. Together, these findings suggest that glial replacement in HD leads to a partial normalization of neuronal gene expression, which is associated with both a restoration of dendritic structural complexity, and a delay in disease progression.

## **Scripts**
Follow the scripts below in order to generate results of this manuscript. 
1. [Preprocessing of individual samples]( https://rawcdn.githack.com/HuynhNPT/hGPC_MSN_rescue/3612e6cb1ee4ef92c0921dfb40b6daae0c7256f7/00_preprocessing_individual_sample.html)</br>

2. [Perform integration]( https://rawcdn.githack.com/HuynhNPT/hGPC_MSN_rescue/3612e6cb1ee4ef92c0921dfb40b6daae0c7256f7/01_integration.html)</br>

3. [Cell type annotation]( https://rawcdn.githack.com/HuynhNPT/hGPC_MSN_rescue/3612e6cb1ee4ef92c0921dfb40b6daae0c7256f7/02_CellTyper.html)</br>

4. Run pySCENIC to derive the gene regulatory network with scripts `03_a_export_for_SCENIC.R`, [03_b_pySCENIC_part1.html]( https://rawcdn.githack.com/HuynhNPT/hGPC_MSN_rescue/3612e6cb1ee4ef92c0921dfb40b6daae0c7256f7/03_b_pySCENIC_part1.html), [03_c_pySCENIC_part2.html]( https://rawcdn.githack.com/HuynhNPT/hGPC_MSN_rescue/3612e6cb1ee4ef92c0921dfb40b6daae0c7256f7/03_c_pySCENIC_part2.html)</br>

5. [Subset integrated object into an object only containing MSN (medium spiny neurons) and re-integrate]( https://rawcdn.githack.com/HuynhNPT/hGPC_MSN_rescue/3612e6cb1ee4ef92c0921dfb40b6daae0c7256f7/04_a_MSN_integration.html)</br>

6. After a new integrated object on MSNs has been produced, we can add the metadata regarding groups of MSNs behavior in wild type and upon transplantation with script `04_b_MSN_metadata.R`
   
7. [Find differential expression and transcription factors of interest]( https://rawcdn.githack.com/HuynhNPT/hGPC_MSN_rescue/25d24666cba70a74c447e4f1927041cd84dff48f/05_Differential_Expression_Analysis.html) </br>
Accompanying scripts: <br>
```
helper_data_wrangling_functions.R
helper_plot_functions.R
05_getDiff_msn.R
05_egos.R
```

## **Data**
Processed `filtered_feature_bc_matrix` can be downloaded from GEO - Accession number: GSE211986
