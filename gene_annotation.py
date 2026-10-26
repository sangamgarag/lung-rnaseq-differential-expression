
import pandas as pd

results = pd.read_csv("differential_expression_results_filtered.csv")

annotation = pd.read_csv(
    "Human.GRCh38.p13.annot.tsv.gz",
    sep="\t",
    compression="gzip"
)

print("Results columns:")
print(results.columns.tolist())

print("\nAnnotation columns:")
print(annotation.columns.tolist())

print("\nResults shape:", results.shape)
print("Annotation shape:", annotation.shape)


import pandas as pd

# Load filtered DESeq2 results
results = pd.read_csv("differential_expression_results_filtered.csv")

# Load gene annotation
annotation = pd.read_csv(
    "Human.GRCh38.p13.annot.tsv.gz",
    sep="\t",
    compression="gzip"
)

# Keep the annotation columns we need
annotation_subset = annotation[
    ["GeneID", "Symbol", "Description", "GeneType"]
].drop_duplicates(subset="GeneID")

# Match results with annotation using GeneID
annotated_results = results.merge(
    annotation_subset,
    on="GeneID",
    how="left",
    validate="many_to_one"
)

# Check the result
print("Rows before annotation:", len(results))
print("Rows after annotation:", len(annotated_results))
print("Genes with symbols:", annotated_results["Symbol"].notna().sum())
print("Genes without symbols:", annotated_results["Symbol"].isna().sum())

# Save annotated results
annotated_results.to_csv(
    "differential_expression_results_annotated.csv",
    index=False
)

print("\nSaved: differential_expression_results_annotated.csv")
