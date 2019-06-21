import pandas as pd
import glob
import os

quant_dir = "/Users/restifo/NIH/data/kallisto_quant"
abundances = glob.glob(quant_dir + "/**/abundance.h5", recursive=True)
df = pd.DataFrame(columns=["sample", "path"])

for i in range(len(abundances)):
    ab = abundances[i]
    df.loc[i] = [os.path.basename(os.path.dirname(ab)), ab]

df.to_csv(quant_dir + "/metadata.tsv", sep="\t", index=False)
