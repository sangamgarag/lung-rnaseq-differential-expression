import pandas as pd

go = pd.read_csv("GO_biological_process_enrichment_corrected.csv")

significant = go[go["FDR"] < 0.05].copy()

significant.to_csv("significant_GO_terms.csv", index=False)

print("Total GO terms tested:", len(go))
print("Significant GO terms:", len(significant))

print("\nTop 10 significant GO terms:")
print(
    significant[
        ["GO_ID", "GO_Process", "Significant_genes_in_term", "FDR"]
    ].head(10).to_string(index=False)
)

print("\nSaved: significant_GO_terms.csv")