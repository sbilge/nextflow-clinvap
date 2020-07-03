#!/usr/bin/env python3

import db_process.process_pharmacogenomics_effects as pharma
import db_process.process_mechanistic_drugs as mechanistic
import table_content.create_content as content
import table_content.prioritize as prior
import table_content.handling as handle
import db_process.process_driver as driver
from collections import defaultdict
import db_process.helpers as helper
import db_process.query as query
import input_process.snv as snv
import functools
from cyvcf2 import VCF
import pandas as pd
import requests
import argparse
import json
import sys
import os


# Command line parameters:
#
# --input-vcf: input vcf
# --output-json: output file name
# --knowledgebase: database dump
# --genome: genome assembly


input_vcf = sys.argv[1]
output_json = sys.argv[2]
genome = sys.argv[3]
knowledgebase = sys.argv[4]

if len(sys.argv) != 5:
    print("usage: {} <input_vcf> <output_json> <genome_assembly_version> <knowledgebase>".format(
        sys.argv[0]))
    exit(1)

KNOWLEDGEBASE = query.read_knowledgebase(knowledgebase)


############
### SNVs ###
############

# PROCESS VCF

# Read in and process VCF file
mvld = snv.parse_vcf(input_vcf)


# Filter and divide dataframe
mvld_high_moderate = snv.get_high_moderate_effect(mvld)
mvld_modifier = snv.get_modifier_effect(mvld)

# Check if mvld_high_moderate dataframe is empty. If it is empty, use mvld_moderate
if not mvld_high_moderate.empty:
    main_mvld = mvld_high_moderate
else:
    if not mvld_modifier.empty:
        handle.empty_high_impact_mvld()
        main_mvld = mvld_modifier

# HGVSp one letter representation
main_mvld = snv.one_letter_repr(main_mvld)


# list of genes
gene_list = query.get_gene_list(main_mvld)

# fields needed from mlvd:
from_mlvd = main_mvld[[
    "HGNC_ID", "SYMBOL", "one_letter_repr", "Consequence"]]

# QUERY DATABASE

#  Get driver gene annotations
driver_ann_mutation = query.get_driver_annotation(
    main_mvld, KNOWLEDGEBASE, "mutation")
# There is a TODO here
# driver_ann_indel = query.get_driver_annotation(
#     main_mvld, KNOWLEDGEBASE, "indel")


# Query database - Get direct drug targets
direct_variant_targets = query.ann_pharm_variant(
    main_mvld, genome, KNOWLEDGEBASE)


# Filter direct drug targets for variant_class
direct_variant_targets_mut = query.filter_variant_class(
    direct_variant_targets, ["MUT", "MUT;MUT", "null"])


# Query database - Get drugs targeting affected gene
gene_targets = query.ann_pharm_gene(main_mvld, KNOWLEDGEBASE, genome)

# Filter direct drug targets for variant_class
gene_targets_mut = query.filter_variant_class(
    gene_targets, ["MUT", "MUT;MUT", "null"])

# Query database - Get mechanistic drug targets, cancer drugs (later - remove the ones found in direct targets and indirect targets)
mechanistic_drug_targets = query.get_mechanistic_drugs(
    main_mvld, KNOWLEDGEBASE)

if (not mechanistic_drug_targets) and (not driver_ann_mutation) and (not gene_targets_mut) and (not direct_variant_targets_mut):
    handle.empty_database_result()


# Check if mechanistic drug targets are empty or not.
if mechanistic_drug_targets:
    mechanistic_flag = True
else:
    mechanistic_flag = False
    handle.empty_mechanistic()


# Check if driver gene annotation dictionary is empty or not.
if driver_ann_mutation:
    driver_flag = True
else:
    driver_flag = False
    handle.empty_driver_annotation()




# PROCESS DATABASE RESULTS


# Update driver dict
# add tumor string to driver dictionary and remove tumor_type list of dictionaries (same function doing both)
for key in driver_ann_mutation:
    helper.add_tumor_string(key, driver_ann_mutation)



# Update driver dict
# collapse info from sources into one representation


if driver_flag:
    for key, annotation in driver_ann_mutation.items():
        role = set()
        source = []
        source_pmid = []
        reference = []
        reference_source = set()
        tumor = []
        for driver_info in annotation:
            # fn updating global role_set
            helper.populate_set(driver_info, "driver_role", role)
            # fn updating global source_list
            helper.populate_list(driver_info, "source_name", source)
            # fn updating global source_pmid_list
            helper.populate_list(driver_info, "source_pmid", source_pmid)
            # fn updating global reference_id list
            helper.populate_list(driver_info, "reference_id", reference)
            # fn updating global reference_source list
            helper.populate_set(
                driver_info, "reference_source", reference_source)
            # fn updating global tumor_list
            helper.populate_list(driver_info, "tumor_list", tumor)

        driver_role = driver.assign_driver_role(role)
        agg_source = "|".join(source)
        agg_pmid = "|".join(source_pmid)
        agg_ref = "|".join(reference)
        agg_ref_source = "|".join(list(reference_source))
        agg_tumor = "|".join([t for t in tumor if t != ""])

        # replace the value of the dict. collapse list into one dict
        driver_ann_mutation[key] = {"driver_role": driver_role, "source_name": agg_source,
                                    "source_pmid": agg_pmid, "reference_id": agg_ref, "reference_source": agg_ref_source, "tumor_list": agg_tumor}

    #  driver information dictionary to dataframe
    df_driver_ann_mutation = driver.dict_to_dataframe(driver_ann_mutation)

else:
    df_driver_ann_mutation = pd.DataFrame(columns=['hgnc_id', 'driver_role', 'source_name', 'source_pmid',
                                                   'reference_id', 'reference_source', 'tumor_list'])

# Start building up driver content. References will be mapped later.
driver_ann = df_driver_ann_mutation[[
    "hgnc_id", "driver_role", "reference_id", "tumor_list"]]

driver_content = pd.merge(from_mlvd, driver_ann, how="right",
                          left_on="HGNC_ID", right_on="hgnc_id").drop(columns=["HGNC_ID"])

# Get references
ref_driver = content.get_references(df_driver_ann_mutation)

###############################################################################################################

# DIRECT VARIANT TARGETS MUT

if direct_variant_targets_mut:
    empty_variant_targets = False

    # DIRECT PHARMACOGENOMICS

    # seperate non-empty pharmacogenomics therapeutics from variant targets
    direct_pharmacogenomics = pharma.get_subset_dict(
        direct_variant_targets_mut, "pharmacogenomics_therapeutics")

    # merge tumor_type list into a list, remove "evidence_statement" and "rating"
    pharma.apply_add_tumor_string(
        direct_pharmacogenomics, "pharmacogenomics_therapeutics")

    # Direct pharmacogenomics effect dictionary to dataframe
    df_direct_pharmacogenomics = pharma.dict_to_dataframe(
        direct_pharmacogenomics, "pharmacogenomics_therapeutics")

    

    if df_direct_pharmacogenomics.empty:
        df_direct_pharmacogenomics = handle.empty_dataframe_direct_pharmacogenomics()

    print()
    # DIRECT PHARMACOGENOMICS FOR COMBINED VARIANTS

    # seperate non-empty pharmacogenomics combined variant therapeutics from variant targets
    direct_pharmacogenomics_combined = pharma.get_subset_dict(
        direct_variant_targets_mut, "pharmacogenomics_combined_variants_therapeutics")

    # Filter based on the distrupted pair gene
    direct_pharmacogenomics_combined = pharma.filter_gene_pair(
        gene_list, direct_pharmacogenomics_combined, "pharmacogenomics_combined_variants_therapeutics")


    # merge tumor_type list into a list, remove "evidence_statement" and "rating"
    pharma.apply_add_tumor_string(
        direct_pharmacogenomics_combined, "pharmacogenomics_combined_variants_therapeutics")

    # Direct pharmacogenomics combined variant effect dictionary to dataframe
    df_direct_pharmacogenomics_combined = pharma.dict_to_dataframe(
        direct_pharmacogenomics_combined, "pharmacogenomics_combined_variants_therapeutics")


    if df_direct_pharmacogenomics_combined.empty:
        pass

    # DIRECT ADVERSE EFFECTS

    # seperate non-empty adverse effect therapeutics from variant targets
    direct_adverse = pharma.get_subset_dict(
        direct_variant_targets_mut, "adverse_effect_therapeutics")

    # Direct adverse effect dictionary to dataframe
    df_direct_adverse = pharma.dict_to_dataframe(
        direct_adverse, "adverse_effect_therapeutics")

    if df_direct_adverse.empty:
        df_direct_adverse = handle.empty_dataframe_direct_adverse()
else:
    handle.empty_direct_variant_annotation()
    empty_variant_targets = True


#############################################################################################################

# GENE TARGETS MUT

if gene_targets_mut:
    empty_gene_targets = False

    # PHARMACOGENOMICS

    # seperate non-empty pharmacogenomics therapeutics from gene targets
    pharmacogenomics_therapeutics = pharma.get_subset_dict(
        gene_targets_mut, "pharmacogenomics_therapeutics")


    # merge tumor_type list into a list, remove "evidence_statement" and "rating"
    pharma.apply_add_tumor_string(
        pharmacogenomics_therapeutics, "pharmacogenomics_therapeutics")

    # Pharmacogenomics effects to dataframe
    df_pharmacogenomics = pharma.dict_to_dataframe(
        pharmacogenomics_therapeutics, "pharmacogenomics_therapeutics")

    if df_pharmacogenomics.empty:
        df_pharmacogenomics = handle.empty_dataframe_pharmacogenomics()

    # PHARMACOGENOMICS COMBINED VARIANTS

    # seperate non-empty pharmacogenomics combined variant therapeutics from gene targets
    pharmacogenomics_combined_variants_therapeutics = pharma.get_subset_dict(
        gene_targets_mut, "pharmacogenomics_combined_variants_therapeutics")


    # Filter based on the distrupted pair gene
    pharmacogenomics_combined_variants_therapeutics = pharma.filter_gene_pair(
        gene_list, pharmacogenomics_combined_variants_therapeutics, "pharmacogenomics_combined_variants_therapeutics")


    # merge tumor_type list into a list, remove "evidence_statement" and "rating"
    pharma.apply_add_tumor_string(pharmacogenomics_combined_variants_therapeutics,
                                "pharmacogenomics_combined_variants_therapeutics")


    df_pharmacogenomics_combined = pharma.dict_to_dataframe(
        pharmacogenomics_combined_variants_therapeutics, "pharmacogenomics_combined_variants_therapeutics")

    if df_pharmacogenomics_combined.empty:
        df_pharmacogenomics_combined = handle.empty_dataframe_pharmacogenomics_combined()


    # ADVERSE EFFECTS

    # seperate non-empty adverse effect therapeutics from gene targets
    adverse_effect_therapeutics = pharma.get_subset_dict(
        gene_targets_mut, "adverse_effect_therapeutics")


    df_adverse = pharma.dict_to_dataframe(
        adverse_effect_therapeutics, "adverse_effect_therapeutics")


    if df_adverse.empty:
        df_adverse = handle.empty_dataframe_adverse()

    # ADVERSE EFFECTS COMBINED VARIANTS

    # seperate non-empty adverse effect combined variants from gene targets
    adverse_effect_combined_variants_therapeutics = pharma.get_subset_dict(
        gene_targets_mut, "adverse_effect_combined_variants_therapeutics")

    # filter based on the distrupted pair
    adverse_effect_combined_variants_therapeutics = pharma.filter_gene_pair(
        gene_list, adverse_effect_combined_variants_therapeutics, "adverse_effect_combined_variants_therapeutics")

    # Adverse effcts combination variants dictionary to dataframe
    df_adverse_combined = pharma.dict_to_dataframe(
        adverse_effect_combined_variants_therapeutics, "adverse_effect_combined_variants_therapeutics")

    if df_adverse_combined.empty:
        df_adverse_combined = handle.empty_dataframe_adverse_combined()
else:
    handle.empty_variant_annotation()
    empty_gene_targets = True



if empty_gene_targets and empty_variant_targets:
    skip_pharmacodynamics_content = True
else:
    skip_pharmacodynamics_content = False


########################################################################################################
if not skip_pharmacodynamics_content:

    # CREATE PHARMACOGENOMICS TABLES CONTENT

    COLUMNS = ["SYMBOL", "drug_name", "evidence_level", "reference_id", "reference_source",
            "variant_drug_association", "tumor_list", "hgnc_id", "variant", "variant_type", "info"]
    COMB_COLUMNS = ["SYMBOL", "drug_name", "evidence_level", "reference_id", "reference_source", "variant_drug_association",
                    "tumor_list", "hgnc_id", "variant", "variant_type", "info", "info_pair"]

    # Content from direct pharmacogenomics effect
    direct_pharm_content = content.get_content(
        df_direct_pharmacogenomics, from_mlvd, COLUMNS, direct=1)


    # Remove direct results covered by df_pharmacogenomics
    pharm = pd.concat(
        [df_pharmacogenomics, df_direct_pharmacogenomics]).drop_duplicates(keep=False)


    # Content from pharmacogenomics effect
    pharm_content = content.get_content(pharm, from_mlvd, COLUMNS)

    # Content from direct pharmacogenomics variant combination effect is expected to be empty for MUT;MUT. No need
    # for trying to integrate this part for [MUT, MUT;MUT]. So skipping now. Will be implemented for CNA and Fusion
    # Content from pharmacogenomics variant combination effect
    pharm_comb_content = content.get_content(
        df_pharmacogenomics_combined, from_mlvd, COMB_COLUMNS, combination=1)

    # move same variant for both genes, if there are any, content to direct pharm table.
    cond = (pharm_comb_content["info"] == "Same gene, same variant, same consequence.") & (
        pharm_comb_content["info_pair"] == "Same gene, same variant, same consequence.")
    rows = pharm_comb_content.loc[cond, :]
    if not rows.empty:
        direct_pharm_content = direct_pharm_content.append(
            rows, ignore_index=True)
        pharm_comb_content.drop(rows.index, inplace=True)


    # Add match level to see whether the results are same gene, same, variant,same consequence or whether the
    # variant and consequence differs in the inut than the database result. For detais see documentation

    sorted_direct_pharm_content = prior.compatibility_sort(
        direct_pharm_content)

    # Join references of the same rows (collapsing same rows with different reference ids)
    cols = ["SYMBOL", "drug_name", "variant_drug_association",
            "tumor_list", "hgnc_id", "variant", "variant_type", "match_level"]
    sorted_direct_pharm_content = prior.unnest_ref(sorted_direct_pharm_content, cols)


    pharm_content = pd.concat([pharm_content, pharm_comb_content])
    sorted_pharm_content = prior.compatibility_sort(pharm_content)
    sorted_pharm_content = prior.unnest_ref(sorted_pharm_content, cols)


    # Get references
    ref_direct_pharm = content.get_references(direct_pharm_content)

    ref_pharm = content.get_references(pharm_content)



    ##############################################################################################


    # CREATE ADVERSE EFFECT TABLE CONTENTS

    ADV_COLUMNS = ["SYMBOL", "drug_name", "evidence_level", "reference_id", "reference_source",
                    "variant_drug_association", "hgnc_id", "variant", "variant_type", "info"]

    ADV_COMB_COLUMNS = ["SYMBOL", "drug_name", "evidence_level", "reference_id", "reference_source", "variant_drug_association",
                        "hgnc_id", "variant", "variant_type", "info", "info_pair"]

    # Content from direct adverse effects table
    direct_adverse_content = content.get_content(
        df_direct_adverse, from_mlvd, ADV_COLUMNS, direct=1, adverse=1)
    

    #remove direct results covered by df_adverse
    adverse = pd.concat([df_direct_adverse, df_adverse]
                        ).drop_duplicates(keep=False)

    # Content from adverse effects table
    adverse_content = content.get_content(
        adverse, from_mlvd, ADV_COLUMNS, adverse=1)


    # Content from adverse effects combined variants
    adverse_comb_content = content.get_content(
        df_adverse_combined, from_mlvd, ADV_COMB_COLUMNS, combination=1, adverse=1)


    # merge direct_adverse_content, adverse_content, adverse_comb_content
    adverse_content = pd.concat(
        [direct_adverse_content, adverse_content, adverse_comb_content]).drop_duplicates()



    # add merge level and sort based on levels
    cols = ["SYMBOL", "drug_name", "variant_drug_association",
            "hgnc_id", "variant", "variant_type", "match_level"]
    sorted_adverse_content = prior.compatibility_sort(adverse_content)
    sorted_adverse_content = prior.unnest_ref(sorted_adverse_content, cols)
    
   
    ref_adverse = content.get_references(adverse_content)

else:
    sorted_direct_pharm_content = handle.empty_pharmacodynamics_place_holder()
    sorted_pharm_content = handle.empty_pharmacodynamics_place_holder()
    sorted_adverse_content = handle.empty_pharmacodynamics_place_holder(adverse=1)
    ref_pharm = handle.empty_reference_place_holder()
    ref_direct_pharm = handle.empty_reference_place_holder()
    ref_adverse = handle.empty_reference_place_holder()
    
    

#############################################################################################################


if mechanistic_flag:

    # Update mechanistic drug targets dict  
    mechanistic_drug_targets = mechanistic.reshape(mechanistic_drug_targets, {})

    # collapse info from sources into one representation
    temp = {}
    for key, annotation in mechanistic_drug_targets.items():
        temp[key] = []
        for db_id, info in annotation.items():
            if len(info) == 1:
                temp[key].append(info[0])
            elif db_id != "null":
                approval = set()
                drug = set()
                source = []
                source_pmid = []
                reference = []
                reference_source = set()
                for drug_info in info:
                    helper.populate_set(drug_info, "approval_status", approval)
                    helper.populate_set(drug_info, "drug_name", drug)
                    helper.populate_list(drug_info, "source_name", source)
                    helper.populate_list(drug_info, "source_pmid", source_pmid)
                    helper.populate_list(drug_info, "reference_id", reference)
                    helper.populate_set(drug_info, "reference_source", reference_source)

                approval_status = mechanistic.assign_approval_status(approval)
                name = list(drug)[0]
                agg_source = "|".join(source)
                agg_pmid = "|".join(source_pmid)
                agg_ref = "|".join(reference)
                agg_ref_source = "|".join(list(reference_source))

                temp[key].append(
                    {"drugbank_id": db_id, "drug_name": name, "approval_status": approval_status,
                        "source_name": agg_source, "reference_id": agg_ref, "reference_source": agg_ref_source,
                        "source_pmid": agg_pmid})
            else:
                for drug_info in info:
                    temp[key].append(drug_info)
    mechanistic_drug_targets = temp

    # filter not-approved drugs from mechanistic drugs
    mechanistic_drugs_filtered = mechanistic.filter_approval_info(
        mechanistic_drug_targets)

    # seperate investigational and approved mechanistic drugs
    approved_mechanistic_drugs = mechanistic.divide_approved_investigational(
        mechanistic_drugs_filtered)[0]


    investigational_mechanistic_drugs = mechanistic.divide_approved_investigational(
        mechanistic_drugs_filtered)[1]


    # Mechanistic drug targets information to dataframe
    if approved_mechanistic_drugs:
        main_mechanistic = approved_mechanistic_drugs
    else:
        handle.empty_approved_mechanistic()
        main_mechanistic = investigational_mechanistic_drugs

    df_mechanistic_targets = mechanistic.dict_to_dataframe(main_mechanistic)
else:
    df_mechanistic_targets = pd.DataFrame(
        columns=['drugbank_id', 'drug_name', 'approval_status', 'source_name', 'reference_id', 'reference_source', 'source_pmid', 'hgnc_id'])


# CREATE MECHANISTIC DRUG TARGETS TABLE CONTENT


if not skip_pharmacodynamics_content:
    # remove the content found in pharmacogenomics therapeutics - bunlardan biri bossa bura patlar mi?
    all_pharm_drugs = pd.concat([direct_pharm_content[["drug_name", "hgnc_id"]], pharm_content[[
        "drug_name", "hgnc_id"]], pharm_comb_content[["drug_name", "hgnc_id"]]]).drop_duplicates()

    cond = df_mechanistic_targets["drug_name"].isin(
        all_pharm_drugs["drug_name"]) & df_mechanistic_targets["hgnc_id"].isin(all_pharm_drugs["hgnc_id"])
    mechanistic = df_mechanistic_targets.loc[~cond]
else:
    mechanistic = df_mechanistic_targets

MECH_COLS = ["SYMBOL", "drug_name", "approval_status",
             "reference_id", "reference_source", "hgnc_id"]
mechanistic_content = content.get_mechanistic_content(
    mechanistic, from_mlvd, MECH_COLS)


ref_mechanistic = content.get_references(mechanistic_content)


##############################################################################################

# CREATE APPENDIX TABLES

# appendix table for all variants # TODO burda bir suphe var
appendix_content = content.appendix_variants(main_mvld)


# Appendix table for references


# concatanate reference dataframes
dataframes = [ref_driver, ref_direct_pharm, ref_pharm, ref_mechanistic, ref_adverse]

ref_df = functools.reduce(
    lambda left, right: pd.concat([left, right], ignore_index=True), dataframes) # concatanate dataframes to create one reference dataframe


# create pubmed query string for eutils
pubmed_ids, ref_map_df = content.id_string(ref_df)
URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?retmode=json;db=pubmed;id={}".format(
    pubmed_ids)


# get API results, create dataframe with reference information
ref_info_df = content.get_ref_details(URL)

# merge reference map and reference info dataframes two get a complete reference map.
ref_map_complete = ref_map_df.merge(ref_info_df, how ="left", left_on="reference_id", right_index=True)

reference_content, reference_map_dict = content.reference_table(ref_map_complete)

# Map reference ids to the dataframes and sort based on the number of the references they have

# get reference map ids for driver_content and sort based n number of references
# listing unknown type of drivers under known typed ones.
driver_content = content.reference_map_ids(driver_content, reference_map_dict)
driver_content = prior.driver_type_sort(driver_content)


# get reference map ids for sorted_direct_pharm_content
sorted_direct_pharm_content = content.reference_map_ids(sorted_direct_pharm_content, reference_map_dict)

# get reference map ids for sorted_pharm_content
sorted_pharm_content = content.reference_map_ids(sorted_pharm_content, reference_map_dict)

# get reference map ids for mechanistic_content and sort based n number of references
mechanistic_content = content.reference_map_ids(mechanistic_content, reference_map_dict)
mechanistic_content = prior.reference_sort(mechanistic_content)

# get reference map ids for sorted_adverse_content
sorted_adverse_content = content.reference_map_ids(sorted_adverse_content, reference_map_dict)


##############################################################################################

# COUNT NUMBERS


# Get number of nonsynonymous genes
num_nonsynonymous = content.nonsynonymous_count(main_mvld)

# Get number of oncogenes
num_oncogene = content.oncogenes_count(driver_content)

# Get number of tumor suppressor genes
num_tsg = content.tsg_count(driver_content)

# Get number of genes, both onco and tumor suppressor genes.
num_onco_tsg = content.onco_tsg_count(driver_content)

# Get unknown type
num_unknown = content.unknown_count(driver_content)


##############################################################################################

# CONVERT CONENT INTO JSON


report_json = {}
driver_table_json = driver_content.to_dict("records")
direct_pharm_table = sorted_direct_pharm_content.to_dict("records")
pharm_table = sorted_pharm_content.to_dict("records")
mechanistic_table = mechanistic_content.to_dict("records")
adverse_table = sorted_adverse_content.to_dict("records")
appendix_variant_table = appendix_content.to_dict("records")
appendix_reference_table = reference_content.to_dict("records")
report_json["num_nonsynonymous"] = num_nonsynonymous
report_json["num_oncogene"] = num_oncogene
report_json["num_tsg"] = num_tsg
report_json["num_onco_tsg"] = num_onco_tsg
report_json["num_unknown"] = num_unknown

report_json["driver_table"] = driver_table_json
report_json["direct_pharm_table"] = direct_pharm_table
report_json["pharm_table"] = pharm_table
report_json["mechanistic_drug_table"] = mechanistic_table
report_json["adverse_table"] = adverse_table
report_json["appendix_variant_table"] = appendix_variant_table
report_json["appendix_reference_table"] = appendix_reference_table


with open(output_json, "w") as f:
    json.dump(report_json, f, indent=4)
