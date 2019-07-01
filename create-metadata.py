import pandas as pd
import glob
import os
import h5py
import datetime as dt

# TODO: Make this whole file part of our snakemake pipeline. This current solution hardcodes the study, which is fine
# TODO: for now and the forseeable future.

months = ["january", "february", "march", "april", "may", "june",
          "july", "august", "september", "october", "november", "december"]
columns = ["tumor_num", "trial", "sample", "path", "method", "response", "sex", "age"]

trial_dirs = ["/data/Robinson-SB/16-C-0027", "/data/Robinson-SB/14-C-0022"]

# Output path for the metadata file
metadata_path = "/data/Robinson-SB/metadata-trials16and14.tsv"

df = pd.DataFrame(columns=columns)
cur_idx = 0
for trial_dir in trial_dirs:
    # response: CR, PR
    # non-response: NR, SD, PD
    this_trial = os.path.basename(trial_dir)
    print("Processing trial", this_trial)

    # Samples taken after this cutoff used RNA access while samples taken before use poly-A
    method_cutoff = dt.datetime(2017, 6, 15)

    quant_dir = os.path.join(trial_dir, "kallisto_quant")

    # Get response data in a DataTable
    resp_data_f = os.path.join(trial_dir, "response_data/%s-tumor-response-sex-age.tsv" % this_trial)
    resp_data = pd.read_csv(resp_data_f, delim_whitespace=True)

    abundances = glob.glob(quant_dir + "/**/abundance.h5", recursive=True)

    for ab in abundances:
        ab_name = os.path.basename(os.path.dirname(ab))

        # Get metadata for this tumor
        tumor_metadata = resp_data[resp_data.tumor == int(ab_name[:4])]

        # Sanity check
        try:
            h5py.File(ab, "r")  # open hdf5 file for reading. If it opens without an OSError, we know it's valid
        except OSError:
            print("Warning: Cannot read HDF5 file %s. Skipping" % ab)
            continue

        if len(tumor_metadata.columns) != 4:
            print("Warning: Sample %s has incomplete metadata, ncol=%d. Skipping" % (ab, len(tumor_metadata.columns)))
            continue

        # Important robustness tip: The tumor metadata might have multiple rows for each tumor, each representing
        # a different follow up. We are just going to take the first row for simplicity at this stage. However, we lose
        # important data when we do this, so TODO: incorporate this into analyses in the future.
        if len(tumor_metadata) > 1:
            tumor_metadata = tumor_metadata[:1]  # take just the first row

        did_respond = "R" if tumor_metadata.response.item() in ["CR", "PR"] else "NR"
        age = int(tumor_metadata.age.item())
        sex = tumor_metadata.sex.item()
        if sex not in ["M", "F"]:
            print("Warning: Sample %s has invalid sex, skipping" % ab)  # even though it's 2019, sex must be m or f
            continue

        tumor_num = tumor_metadata.tumor.item()

        # Need to convert this string into a datetime. Problem is, they are formatted differently but all contain a
        # [month]_[day]_[year] string. The following code searches for the beginning of this string by looking
        # for the month string
        date_str = os.path.basename(os.path.dirname(os.path.dirname(ab))).lower()
        dt_obj = None
        for month in months:
            try:
                start = date_str.index(month)
                parsed = "_".join(date_str[start:].split("_")[:3])

                # TODO: If parsing with this format string fails, try another one instead of just skipping the sample
                dt_obj = dt.datetime.strptime(parsed, "%B_%d_%Y")
            except ValueError:
                pass

        if dt_obj is None:
            print("Warning: Couldn't identify date for sample %s, assuming method is poly-A" % ab)

            # by giving a date 1 day before the cut off, this sample's method will be poly-A
            dt_obj = method_cutoff - dt.timedelta(days=1)

        method = "rna_access" if dt_obj > method_cutoff else "polya"
        df.loc[cur_idx] = [tumor_num, this_trial, ab_name, ab, method, did_respond, sex, age]
        cur_idx += 1


print("Writing metadata file")
df.to_csv(metadata_path, sep="\t", index=False)
