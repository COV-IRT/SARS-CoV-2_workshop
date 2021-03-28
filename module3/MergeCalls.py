import argparse
import vcf
import pandas as pd

def read_ivar(filename):
    ivar_calls = pd.read_csv(filename, sep='\t')
    ivar_calls["Variant"] = ivar_calls.apply(lambda row: 
                                     str(row["POS"]) + \
                                     str(row["REF"]) + ">" + \
                                     str(row["ALT"]), axis=1)
    ivar_calls = ivar_calls.drop_duplicates(subset=["Variant"])

    return ivar_calls


def read_lofreq(filename):
    lofreq_calls = pd.DataFrame(columns=["CHROM", "POS", "REF", "ALT", "QUAL", 
                                     "REF_DP", "REF_RV", "ALT_DP", "ALT_RV",
                                     "ALT_FREQ", "TOTAL_DP"])
    vcf_reader = vcf.Reader(filename=filename)
    for row in vcf_reader:
        lofreq_calls = lofreq_calls.append({"CHROM": row.CHROM, 
                                            "POS": int(row.POS),
                                            "REF": row.REF,
                                            "ALT": row.ALT[0],
                                            "QUAL": row.QUAL, 
                                            "REF_DP": row.INFO["DP4"][0] + row.INFO["DP4"][1], 
                                            "REF_RV": row.INFO["DP4"][1],
                                            "ALT_DP": row.INFO["DP4"][2] + row.INFO["DP4"][3],
                                            "ALT_RV": row.INFO["DP4"][3],
                                            "ALT_FREQ": row.INFO["AF"],
                                            "TOTAL_DP": row.INFO["DP"],
                                           }, 
                                           ignore_index=True)

    lofreq_calls["Variant"] = lofreq_calls.apply(lambda row: 
                                         str(row["POS"]) + \
                                         str(row["REF"]) + ">" + \
                                         str(row["ALT"]), axis=1)
    return lofreq_calls


def main():
    parser = argparse.ArgumentParser(description="Variant call merging for LoFreq and iVar")
    parser.add_argument("-m", "--min-af", type=float, default=0.02, 
        help="Minimum allele frequncy of the variants (between 0 and 1)",
        required=True)
    parser.add_argument("-p", "--pass-only", action="store_true",
        help="Only retain vairants with p-value <= 0.05 from iVar",
        required=False)
    parser.add_argument("ivar_input", type=str, action="store",
        help="Path to the iVar output .tsv")
    parser.add_argument("lofreq_input", type=str, action="store",
        help="Path to the LoFreq output .vcf")
    parser.add_argument("-o", "--output", type=str, default="output.tsv",
        help="Name of the output file for merged calls")
    args = parser.parse_args()

    ivar_calls = read_ivar(args.ivar_input)
    lofreq_calls = read_lofreq(args.lofreq_input)
    merged_calls = lofreq_calls.merge(ivar_calls, on=["Variant"], suffixes=("_LoFreq", "_iVar"))
    filtered_merged_calls = merged_calls[~((merged_calls["ALT_FREQ_iVar"] < args.min_af) & \
                                           (merged_calls["ALT_FREQ_LoFreq"] < args.min_af))]
    if args.pass_only:
        filtered_merged_calls = filtered_merged_calls[filtered_merged_calls["PASS"] == True]

    filtered_merged_calls.to_csv(args.output, sep="\t", index=False)

if __name__ == "__main__":
    main()
