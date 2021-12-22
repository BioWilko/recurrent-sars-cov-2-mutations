import pandas as pd
import re
import sys
import argparse
import pathlib

parser = argparse.ArgumentParser()
parser.add_argument("-l", "--long_df", help="Long df CSV", type=str)
parser.add_argument(
    "-cm",
    "--cluster_mutations",
    help="Cluster mutations by AA position",
    action="store_true",
)
parser.add_argument("-p", "--prefix", help="prefix for output file", type=str)
parser.add_argument("-o", "--outdir", help="outpit directory", type=pathlib.Path)
parser.add_argument(
    "-c", "--categorise", help="Categorise patients", action="store_true"
)
parser.add_argument(
    "-e",
    "--export_occurences",
    help="Export the occurence df generated as an intermediate to a CSV",
    action="store_true",
)
args = parser.parse_args()

long_df = pd.read_csv(args.long_df)

spike_locs = pd.read_csv(
    "/home/sam/onedrive/bioinformatics/mutation_rate_metastudy/spike_map.csv"
)

pt_categories = pd.read_csv(
    "/home/sam/onedrive/bioinformatics/mutation_rate_metastudy/dataset/patient_categories.csv",
    index_col="pt",
)

pt_category_list = set(pt_categories["category"]) if args.categorise else ["n/a"]

category_filename_dict = {
    "Combined": "C",
    "T cell dominant": "T",
    "B cell dominant": "B",
}

long_df[["pt", "day"]] = long_df["sample"].str.split("_", 1, expand=True)

long_df["day"] = long_df["day"].str.split("_", expand=True)[0]

long_df["day"] = long_df["day"].astype(int)

if args.categorise:
    for index, row in long_df.iterrows():
        long_df.at[index, "category"] = pt_categories.at[row["pt"], "category"]


desired_cols = (
    [
        "pt",
        "category",
        "day",
        "refpos",
        "varclass",
        "refvar",
        "qvar",
        "protein",
        "annotation",
        "varname",
    ]
    if args.categorise
    else [
        "pt",
        "day",
        "refpos",
        "varclass",
        "refvar",
        "qvar",
        "protein",
        "annotation",
        "varname",
    ]
)

long_df = long_df[desired_cols]


pt_groups = long_df.groupby("pt")

occurence_df = pd.DataFrame()

col_names = (
    [
        "pt",
        "category",
        "day",
        "ref_pos",
        "mut_type",
        "ref_nt",
        "alt_nt",
        "protein",
        "protein_annotation",
        "annotation",
        "facet_annotation",
    ]
    if args.categorise
    else [
        "pt",
        "day",
        "ref_pos",
        "mut_type",
        "ref_nt",
        "alt_nt",
        "protein",
        "protein_annotation",
        "annotation",
        "facet_annotation",
    ]
)

for pt in pt_groups.groups:
    df = pt_groups.get_group(pt).copy()
    df.sort_values(by="day", inplace=True)
    df = df.query(
        "varclass != 'extragenic' & varclass != 'SNP_silent' & varclass != 'deletion_frameshift'"
    )
    called_mutations = []
    day_0_mutations = df.loc[df["day"] == 0]["varname"].to_list()
    pt_df_lists = []
    for index, row in df.iterrows():
        if row["day"] != 0:
            if (
                row["varname"] not in day_0_mutations
                and row["varname"] not in called_mutations
            ):
                append_list = row.to_list()
                if str(row["varname"]).startswith("S:"):
                    for idx, r in spike_locs.iterrows():
                        codon = int(re.sub("[^0-9]", "", row["varname"]))
                        if codon >= r["start"] and codon <= r["end"]:
                            append_list.append(r["domain"])
                            break
                        elif idx >= 9:
                            append_list.append("other")
                elif str(row["varname"]).split(":")[0].startswith("NSP"):
                    append_list.append("ORF1ab")
                else:
                    append_list.append("not_s/o")
                pt_df_lists.append(append_list)
                called_mutations.append(row["varname"])
    # print(pt_df_lists)
    append_df = pd.DataFrame(pt_df_lists, columns=col_names)
    occurence_df = occurence_df.append(append_df, ignore_index=True)

if args.export_occurences:
    occurence_df.to_csv(f"{args.outdir}/{args.prefix}_occurences.csv", index=False)

deletion_df = occurence_df.loc[occurence_df["mut_type"] == "deletion"]
deletion_df = deletion_df.sort_values(by="ref_pos").copy()
deletions_clustered = []
for index, row in deletion_df.iterrows():
    if row["annotation"] not in deletions_clustered:
        mask = (deletion_df["ref_pos"] >= row["ref_pos"]) & (
            deletion_df["ref_pos"] <= row["ref_pos"] + 18
        )
        cluster = deletion_df.loc[mask]
        annotation_split = row["annotation"].split(":")
        codon = int(re.sub("[^0-9]", "", annotation_split[1]))
        delta_annotation = f"{annotation_split[0]}:Δ{codon} region"
        for idx, r in cluster.iterrows():
            occurence_df.at[idx, "annotation"] = delta_annotation
            deletion_df.drop(idx, inplace=True)
            deletions_clustered.append(delta_annotation)


def cluster_annotations(annotation):
    if annotation.endswith("region"):
        return annotation
    else:
        return annotation[:-1] + "X"


if args.cluster_mutations:
    clustered_annotations = occurence_df["annotation"].apply(cluster_annotations)
    occurence_df["annotation"] = clustered_annotations

occurence_df["day"] = occurence_df["day"].astype(int)


def cumulative_count(dataframe, mutation_list):
    cumulative_rows = []
    for day in range(0, dataframe["day"].max() + 1):
        to_date = dataframe.loc[dataframe["day"] <= day]
        if len(to_date) != 0:
            for index, row in mutation_list.iterrows():
                count = len(
                    to_date.loc[to_date["annotation"] == str(row["annotation"])]
                )
                append_list = [day]
                append_list.extend(row.to_list())
                append_list.append(count)
                cumulative_rows.append(append_list)
        else:
            for index, row in mutation_list.iterrows():
                count = 0
                append_list = [day]
                append_list.extend(row.to_list())
                append_list.append(count)
                cumulative_rows.append(append_list)

    cumulative_df = pd.DataFrame(
        cumulative_rows,
        columns=[
            "day",
            "ref_pos",
            "mut_type",
            "ref_nt",
            "alt_nt",
            "protein",
            "protein_annotation",
            "annotation",
            "facet_annotation",
            "count",
        ],
    )

    cumulative_df.drop_duplicates(
        subset=["annotation", "day"], keep="first", inplace=True
    )
    return cumulative_df


if args.categorise:
    for category in pt_category_list:

        category_df = occurence_df.loc[occurence_df["category"] == category]

        mutations = category_df[
            [
                "ref_pos",
                "mut_type",
                "ref_nt",
                "alt_nt",
                "protein",
                "protein_annotation",
                "annotation",
                "facet_annotation",
            ]
        ].copy()

        mutation_list = mutations.drop_duplicates(ignore_index=True)

        category_cumulative = cumulative_count(category_df, mutation_list)

        category_cumulative.to_csv(
            f"{args.outdir}/{args.prefix}_{category}.csv", index=False,
        )
else:
    mutations = occurence_df[
        [
            "ref_pos",
            "mut_type",
            "ref_nt",
            "alt_nt",
            "protein",
            "protein_annotation",
            "annotation",
            "facet_annotation",
        ]
    ].copy()

    mutation_list = mutations.drop_duplicates(ignore_index=True)

    cumulative = cumulative_count(occurence_df, mutation_list)

    cumulative.to_csv(f"{args.outdir}/{args.prefix}_cumulative.csv", index=False)
