#!/usr/bin/env python3

import json
import sys
import pandas as pd
import copy
import os


main_json = sys.argv[1]
metadata_json = sys.argv[2]
output_json = sys.argv[3]
do = sys.argv[4]



if len(sys.argv) != 5:
    print("usage: {} <main_json> <metadata_json> <output_json> <disease_ontology_map>".format(
        sys.argv[0]))
    exit(1)


# Read metadata file
with open(metadata_json) as m:
    metadata = json.load(m)

# Read main report file
with open(main_json) as r:
    report = json.load(r)

# merge main report content with metadata
merge_report = metadata.copy()
merge_report.update(report)

# Load disease vocabulary
dtypes = {"tumor_type": "str", "disease_name": "str",
          "doid": "str", "icd10": "str", "do_name": "str", "abbr": "str"}
vocab = pd.read_csv(do, sep="\t", dtype=dtypes)
vocab = vocab.apply(lambda x: x.str.lower())


def get_field(key, metadata):
    val = metadata.get(key, "null")
    if not val:
        val = "null"
    return val


def give_na_diagnosis_tag(report, *keys):
    copy_report = copy.deepcopy(report)
    for key in keys:
        for d in copy_report[key]:
            d["seen_in_diagnosis"] = "Not applicable"
    return copy_report


diagnosis = get_field("diagnosis", merge_report).lower()
do_name = get_field("do_name", merge_report).lower()
doid = get_field("doid", merge_report).lower()
icd10 = get_field("icd10", merge_report).lower()
abbr = get_field("abbreviation", merge_report).lower()


if all(val == "null" for val in [diagnosis, do_name, doid, icd10, abbr]):
    tagged_report = give_na_diagnosis_tag(merge_report, "driver_table",
                                          "direct_pharm_table", "pharm_table")
    with open(output_json, "w") as o:
        json.dump(tagged_report, o, indent=4)
    sys.exit(0)


def in_vocabulary(vocab_str, column_name, result_list):
    vocab_list = vocab_str.split(",")
    for v in vocab_list:
        if vocab_str == "icd10":
            if not vocab[vocab[column_name].str.contains(v)].empty:
                result_list.append(v)
    else:
        if not vocab[vocab[column_name] == v].empty:
            result_list.append(v)
    return result_list


def get_tumor_name(vocab_list, vocab_df, vocab_df_column_name):
    """Function to obtain tumor type keywords as given databases from disease ontology id, icd10 codes or do_"""
    conv_list = []
    for v in vocab_list:
        cond = vocab_df[vocab_df_column_name] == v
        indices = cond[cond].index
        conv_list = vocab_df.loc[indices, "tumor_type"].tolist()
    return conv_list


vocab_diagnosis = []
vocab_diagnosis2 = []
vocab_diagnosis3 = []
vocab_doid = []
vocab_icd10 = []
vocab_abbr = []


if diagnosis != "null":
    vocab_diagnosis = in_vocabulary(diagnosis, "tumor_type", vocab_diagnosis)
    vocab_diagnosis2 = in_vocabulary(
        diagnosis, "disease_name", vocab_diagnosis2)
    vocab_diagnosis3 = in_vocabulary(diagnosis, "do_name", vocab_diagnosis3)
    vocab_diagnosis2_tn = get_tumor_name(
        vocab_diagnosis2, vocab, "disease_name")
    vocab_diagnosis3_tn = get_tumor_name(vocab_diagnosis3, vocab, "do_name")
else:
    vocab_diagnosis2_tn = []
    vocab_diagnosis3_tn = []


if doid != "null":
    vocab_doid = in_vocabulary(doid, "doid", vocab_doid)
    vocab_doid_tn = get_tumor_name(vocab_doid, vocab, "doid")
else:
    vocab_doid_tn = []


if icd10 != "null":
    vocab_icd10 = in_vocabulary(icd10, "icd10", vocab_icd10)
    vocab_icd10_tn = get_tumor_name(vocab_icd10, vocab, "icd10")
else:
    vocab_icd10_tn = []


if abbr != "null":
    vocab_abbr = in_vocabulary(abbr, "abbr", vocab_abbr)
    vocab_abbr_tn = get_tumor_name(vocab_abbr, vocab, "abbr")
else:
    vocab_abbr_tn = []


keyword_lists = [vocab_diagnosis, vocab_diagnosis2_tn,
                 vocab_diagnosis3_tn, vocab_doid_tn, vocab_icd10_tn, vocab_abbr_tn]
if all(result == [] for result in keyword_lists):
    tagged_report = give_na_diagnosis_tag(
        merge_report, "driver_table", "direct_pharm_table", "pharm_table")
    with open(output_json, "w") as o:
        json.dump(tagged_report, o, indent=4)
    sys.exit(0)
else:
    keywords = []
    for x in keyword_lists:
        keywords.extend(x)
    keywords = list(set(keywords))


def tumor_type_match(section, keywords):
    for d in section:
        d["seen_in_diagnosis"] = False
        for k in keywords:
            if k in d["db_tumor_repr"].lower():
                d["seen_in_diagnosis"] = True
                break
            else:
                continue
    return section


merge_report["driver_table"] = tumor_type_match(
    merge_report["driver_table"], keywords)
merge_report["direct_pharm_table"] = tumor_type_match(
    merge_report["direct_pharm_table"], keywords)
merge_report["pharm_table"] = tumor_type_match(
    merge_report["pharm_table"], keywords)


with open(output_json, "w") as o:
    json.dump(merge_report, o, indent=4)
