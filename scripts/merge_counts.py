import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("--samples", dest="samplesheet", default="config/samples.csv")
parser.add_argument("--out", dest="out", default="results/counts/combined_counts.tsv")
args, _ = parser.parse_known_args()

samples = pd.read_csv(args.samplesheet)
mat = None
for _, row in samples.iterrows():
    sample = row["sample"]
    path = f"results/counts/{sample}.counts.txt"
    df = pd.read_csv(path, sep="\t", comment="#")
    # featureCounts format: first 6 columns are gene meta; last column is counts
    df = df[["Geneid", df.columns[-1]]].copy()
    df.columns = ["Geneid", sample]
    if mat is None:
        mat = df
    else:
        mat = mat.merge(df, on="Geneid", how="outer")

mat = mat.fillna(0)
mat.to_csv(args.out, sep="\t", index=False)