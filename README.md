# Lung RNA-seq Differential Expression Analysis

**Dataset:** GSE229705  
**Tools:** R, DESeq2, Python, Pandas, SciPy, Statsmodels

## Project Overview
This project investigates gene-expression differences between lung tumour tissue and adjacent-normal lung tissue using publicly available bulk RNA-seq data.

## Dataset
- 244 samples: 122 tumour and 122 adjacent-normal samples
- 122 matched patient pairs
- 25,002 genes retained after low-expression filtering

## Methodology
1. Prepared and matched gene-count data with sample metadata.
2. Filtered low-expression genes.
3. Performed paired differential-expression analysis using DESeq2 with the design `~ Patient + Tissue`.
4. Annotated results with human gene information.
5. Generated PCA, volcano and heatmap visualizations.
6. Performed Gene Ontology Biological Process enrichment using a hypergeometric test and false-discovery-rate correction.

## Key Results
- 25,002 genes analyzed.
- 17,850 genes had adjusted p-values below 0.05.
- 5,107 genes had adjusted p-values below 0.05 and absolute log2 fold change of at least 1.
- 12 GO Biological Process terms had FDR below 0.05.

## Visualizations
- PCA plot: `PCA_plot_filtered.png`
- Volcano plot: `volcano_plot_filtered.png`
- Top-30-gene heatmap: `heatmap_top30_genes_filtered.png`

## Interpretation
The enriched GO terms included immune-related and translation-related biological processes. Because the dataset represents bulk tissue, these findings may reflect contributions from multiple cell types and should not be interpreted as proof of pathway activation in tumour cells.

## Limitations
Statistical significance does not necessarily imply a large biological effect. Findings require interpretation alongside fold changes, sample-level patterns and the original study context.

## Skills Demonstrated
RNA-seq analysis, paired experimental design, differential-expression analysis, gene annotation, statistical testing, multiple-testing correction, data visualization and biological interpretation.
