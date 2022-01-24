# Supplementary material for the manuscript 'Pneumococcal within-host diversity during colonisation, transmission and treatment'

### Data

- Due to its size, the data folder is attached to the release of this repository.
- The accession numbers to access the raw sequencing data are given in Supplementary Table 1 which is available in the data folder and as an attachment with the manuscript.

### Code

The following code is made available to aid with the reproducibility of our analyses. Due to the large computational requirements and the need to protect the privacy of the participants only some of this code is set up to run out of the box. If any assistance is required to reproduce parts of our analyses please send me an [email](g.tonkinhill@gmail.com)

* **initial_qc.Rmd**
  * The initial quality control filters used to determine sensible thresholds for the inclusion of samples in downstream analyses.
* **lineage_deconvolution_serotype_calling.Rmd**
  * The code used to run and analyse the mGEMS deconvolution pipeline and serotype calling pipelines.
* **variant_calling.Rmd**
  * Code used to perform the within-host SNP variant calling on samples involving only a single lineage.
* **replicate_verification.Rmd**
  * Code used to compare the results of the deconvolution pipeline and variant calling pipelines on the set of samples for which replicate culture and sequencing steps were performed.
* **resistance_and_virulence_calls.Rmd**
  * Code used to call and anlayse resistance and virulence elements.
* **antimicrobial_treatment_pairwise_gwas.Rmd**
  * Code used to perform the antimicrobial treatment GWAS using paired samples taken from the same person before and after treatment.
* **antimicrobial_treatment_pyseer_gwas.Rmd**
  * Code to run the Pyseer GWAS pipeline on samples that had and had not recently been treated with antibiotics.
* **clearance_treatment.Rmd**
  * Analysis code to consider the impacts of antimicrobial treatment on the clearance of lineages and in particular the clearance of PBP gene types.
* **mutation_and_selection.Rmd**
  * Analysis code to characterise selection at the genome and gene level and to additionally compare the mutational spectrum between mutations identified within a single host to that seen across larger evolutionary timescales.
* **pairwise_transmission.Rmd**
  * Code to run the fastTransCluster algorithm and process the resulting transmission estimates.
  
### Scripts

This folder includes various custom scripts that are called by the Rmarkdown files described above.


