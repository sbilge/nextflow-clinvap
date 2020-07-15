#!/usr/bin/env python3
import pandas as pd
import copy
import json
import sys
import os

main_json = sys.argv[1]
output_json = sys.argv[2]
do = sys.argv[3]


if len(sys.argv) != 4:
    print("usage: {} <main_json> <output_json> <disease_ontology_map>".format(
        sys.argv[0]))
    exit(1)

# Load disease vocabulary
dtypes = {"tumor_type": "str", "disease_name": "str", "doid": "str", "icd10": "str", "do_name": "str", "abbr": "str"}
vocab = pd.read_csv(do, sep = "\t", dtype = dtypes)
vocab = vocab.apply(lambda x: x.str.lower())

# Read main report file
with open(main_json) as r:
    report = json.load(r)


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


diagnosis = get_field("diagnosis", report).lower()
do_name = get_field("do_name", report).lower()
doid = get_field("doid", report).lower()
icd10 = get_field("icd10", report).lower()
abbr = get_field("abbreviation", report).lower()


if all(val == "null" for val in [diagnosis, do_name, doid, icd10, abbr]):
    tagged_report = give_na_diagnosis_tag(report, "driver_table",
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
        cond = vocab_df[vocab_df_column_name]==v
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
    vocab_diagnosis2 = in_vocabulary(diagnosis, "disease_name", vocab_diagnosis2)
    vocab_diagnosis3 = in_vocabulary(diagnosis, "do_name", vocab_diagnosis3)
    vocab_diagnosis2_tn = get_tumor_name(vocab_diagnosis2, vocab, "disease_name")
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
    vocab_icd10_tn=[]


if abbr != "null":
    vocab_abbr = in_vocabulary(abbr, "abbr", vocab_abbr)
    vocab_abbr_tn = get_tumor_name(vocab_abbr, vocab, "abbr")
else:
    vocab_abbr_tn = []


keyword_lists = [vocab_diagnosis, vocab_diagnosis2_tn, vocab_diagnosis3_tn, vocab_doid_tn, vocab_icd10_tn, vocab_abbr_tn]
if all(result == [] for result in keyword_lists):
    tagged_report = give_na_diagnosis_tag(report, "driver_table", "direct_pharm_table", "pharm_table")
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


report["driver_table"] = tumor_type_match(report["driver_table"], keywords)
report["direct_pharm_table"] = tumor_type_match(report["direct_pharm_table"], keywords)
report["pharm_table"] = tumor_type_match(report["pharm_table"], keywords)


# Filter report based on diagnosis

def filter_on_diagnosis(section):
    filtered_section = [
        d for d in section if d["seen_in_diagnosis"] != False]
    return filtered_section


report["driver_table"] = filter_on_diagnosis(report["driver_table"])
report["direct_pharm_table"] = filter_on_diagnosis(report["direct_pharm_table"])
report["pharm_table"] = filter_on_diagnosis(report["pharm_table"])



with open(output_json, "w") as o:
    json.dump(report, o, indent=4)


