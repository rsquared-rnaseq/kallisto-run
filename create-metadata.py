import pandas as pd
import glob
import os
from datetime import datetime

months = ["january", "february", "march", "april", "may", "june",
          "july", "august", "september", "october", "november", "december"]

# Samples taken after this cutoff used RNA access while samples taken before use poly-A
method_cutoff = datetime(2017, 6, 15)

quant_dir = "/Users/restifo/NIH/data/kallisto_quant"
abundances = glob.glob(quant_dir + "/**/abundance.h5", recursive=True)
df = pd.DataFrame(columns=["sample", "path", "method"])

for i in range(len(abundances)):
    ab = abundances[i]

    # Need to convert this string into a datetime. Problem is, they are formatted differently but all contain a
    # [month]_[day]_[year] string. The following code searches for the beginning of this string by looking
    # for the month string
    date_str = os.path.basename(os.path.dirname(os.path.dirname(ab))).lower()
    dt_obj = None
    for month in months:
        try:
            start = date_str.index(month)
            parsed = "_".join(date_str[start:].split("_")[:3])
            # print(date_str, " ", start, " ", parsed)

            dt_obj = datetime.strptime(parsed, "%B_%d_%Y")
        except ValueError:
            pass

    if dt_obj is None:
        raise ValueError("Couldn't identify date for sample %s" % ab)

    # TODO: >= might not be correct. get clarification on whether that day is included and if i need to worry about the
    # time of day
    method = "rna_access" if dt_obj >= method_cutoff else "polya"
    df.loc[i] = [os.path.basename(os.path.dirname(ab)), ab, method]

df.to_csv(quant_dir + "/metadata.tsv", sep="\t", index=False)
