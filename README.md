# **Variant-Specific Landscape of Mutual Exclusivity Among BRAF, EGFR, and KRAS Oncogenes in Human Cancer**

Welcome to the repository for our study. This study is currently in the publication process and can be viewed in advance at the following link:
https://www.medrxiv.org/content/10.1101/2023.10.21.23297089v1

## **Summary**
This study explores the mutual exclusivity and co-occurrence patterns of BRAF, KRAS, and EGFR mutations in human cancer. Using automated analysis on diverse datasets, we identified specific mutation subtypes with distinct clinical implications, revealing a lower likelihood of co-occurrence among certain mutations. The findings justify the refinement of mutational classifications and provide a variant-specific database for precision oncology, offering insights that may guide the discovery of novel synthetically lethal interactions for targeted cancer therapy.

## **Structure of the repository**
**Raw data cell lines:**

The input data and parameters for analyzing the cell line data are stored here. The following is a brief overview:
- alterations_across_cell lines.tsv: The mutation data of 1570 different cell lines from the Cancer Cell Line Encyclopedia (Broad, 2019) study of cBioPortal (dataset was current as of December 01, 2022).
- KRAS_Mutation_Classes.txt: Class assignment of certain KRAS mutations
- BRAF_Mutation_Classes.txt: Class assignment of certain BRAF mutations
- EGFR_Mutation_Classes.txt: Class assignment of certain EGFR mutations
- EGFR_Mutation_Classes_v2_complete.txt: Corresponds to EGFR_Mutation_Classes.txt, except that all cases of exon19 deletions have been replaced by the specific deletions in the input data
- Color_Code.txt: Color assignment for the final figures

**Raw data samples:**

The input data and parameters for analyzing the patient data are stored here. The following is a brief overview:
- 2022_12_alterations_across_samples.tsv: The mutation data of 68479 different samples of cBioPortal (dataset was current as of November 15, 2022).
- KRAS_Mutation_Classes.txt: Class assignment of certain KRAS mutations
- BRAF_Mutation_Classes.txt: Class assignment of certain BRAF mutations
- EGFR_Mutation_Classes.txt: Class assignment of certain EGFR mutations
- EGFR_Mutation_Classes_v2_complete.txt: Corresponds to EGFR_Mutation_Classes.txt, except that all cases of exon19 deletions have been replaced by the specific deletions in the input data
- Color_Code.txt: Color assignment for the final figures

**Code:**

- Cell line_Data_ProjectCooccurrence.R: This R code is used to automatically analyze cell line mutation data for cooccurrence and mutual exclusivity. The input folder **'Raw data cell lines '** must be selected, the results are saved in the output folder **'results_cell_lines '**.

- PatientData_ProjectCooccurrence.R: This R code is used to automatically analyze cell line mutation data for cooccurrence and mutual exclusivity. The input folder **'Raw data samples '** must be selected, the results are saved in the output folder **'results_samples '**.

**Summary of sample data:**

In the file '**UniquePatientsData.tsv**' the used patient sample data from cBioPortal is summarized as described in the manuscript after filtering for duplicates etc., assignment of classes for specific gene mutations and with the respective cancer entity.
