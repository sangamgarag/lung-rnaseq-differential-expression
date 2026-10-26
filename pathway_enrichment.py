
import pandas as pd
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests

# Load results and annotation
results = pd.read_csv("differential_expression_results_annotated.csv")
annotation = pd.read_csv(
    "Human.GRCh38.p13.annot.tsv.gz",
    sep="\t",
    compression="gzip",
    low_memory=False
)

# Define tested and significant genes as strings
results["GeneID"] = results["GeneID"].astype(str)
annotation["GeneID"] = annotation["GeneID"].astype(str)

tested = set(results["GeneID"].dropna())
significant = set(
    results.loc[
        results["padj"].notna() & (results["padj"] < 0.05),
        "GeneID"
    ].dropna()
)

# Parse GO IDs and descriptions as aligned pairs
rows = []
mismatched_rows = 0

for _, row in annotation[
    ["GeneID", "GOProcessID", "GOProcess"]
].dropna(subset=["GOProcessID", "GOProcess"]).iterrows():

    gene = row["GeneID"]
    if gene not in tested:
        continue

    ids = str(row["GOProcessID"]).split("///")
    names = str(row["GOProcess"]).split("///")

    if len(ids) != len(names):
        mismatched_rows += 1
        continue

    for go_id, go_name in zip(ids, names):
        go_id = go_id.strip()
        go_name = go_name.strip()

        if go_id and go_name:
            rows.append((gene, go_id, go_name))

go = pd.DataFrame(rows, columns=["GeneID", "GO_ID", "GO_Process"])
go = go.drop_duplicates()

# Over-representation analysis using tested genes as background
N = len(tested)
n = len(significant)
output = []

for go_id, group in go.groupby("GO_ID"):
    term_genes = set(group["GeneID"]) & tested
    overlap = term_genes & significant

    K = len(term_genes)
    k = len(overlap)

    if K == 0 or k == 0:
        continue

    output.append({
        "GO_ID": go_id,
        "GO_Process": group["GO_Process"].iloc[0],
        "Genes_in_background": K,
        "Significant_genes_in_term": k,
        "P_value": hypergeom.sf(k - 1, N, K, n),
        "Overlapping_GeneIDs": ";".join(sorted(overlap))
    })

enrichment = pd.DataFrame(output)

if enrichment.empty:
    print("No GO terms with overlapping significant genes were found.")
else:
    enrichment["FDR"] = multipletests(
        enrichment["P_value"], method="fdr_bh"
    )[1]

    enrichment = enrichment.sort_values("FDR")
    enrichment.to_csv("GO_biological_process_enrichment_corrected.csv",
                      index=False)

    print("Tested genes:", N)
    print("Significant genes:", n)
    print("Rows skipped due to ID/description count mismatch:",
          mismatched_rows)
    print("GO terms tested:", len(enrichment))
    print("GO terms with FDR < 0.05:",
          (enrichment["FDR"] < 0.05).sum())
    print("\nTop 15 GO terms:")
    print(enrichment[
        ["GO_ID", "GO_Process", "Significant_genes_in_term", "P_value", "FDR"]
    ].head(15).to_string(index=False))

    print("\nSaved: GO_biological_process_enrichment_corrected.csv")
