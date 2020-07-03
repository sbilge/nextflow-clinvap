### VALID ONLY FOR SNVs

import json
import pandas
import sys
import copy


def read_knowledgebase(filename):
    """Helper fn, Read Knowledgebase"""
    with open(filename) as db:
        cancer_db = json.load(db)
    return cancer_db


def get_gene_list(dataframe):
    """Helper fn, function to get gene list from the dataframe"""
    new_df = dataframe[["HGNC_ID"]]
    new_df = new_df.drop_duplicates()
    gene_list = new_df["HGNC_ID"].tolist()
    return gene_list


def get_inner_dict(key1, key2, dictionary):
    """Helper fn, function to get second level inner dictionary from a nested one"""
    if key1 in dictionary:
        inner_dict = dictionary[key1].get(key2)
        inner_d = copy.deepcopy(inner_dict)
        if inner_dict:
            return inner_d
        else:
            return None


def get_variant_list(dataframe):
    """Helper fn, get variant coordinates into a list (dataframe[fields].drop_dupicates might be necessary)"""
    fields = ["HGNC_ID", "chr", "start", "end", "ref", "alt"]
    rows = dataframe[fields].values.tolist()
    return rows


def get_driver_annotation(dataframe, knowledgebase, v_type):
    """Get driver gene annotation, query with hgnc is and variant class. Returns driver genes"""
    driver_mutation_annotation = {}
    driver_indel_annotation = {}
    gene_list = get_gene_list(dataframe)
    for gene in gene_list:
        driver = get_inner_dict(gene, "driver_annotation", knowledgebase)
        if not isinstance(driver, dict):
            continue
        mutation = driver.get("mutation")
        indel = driver.get("indel")
        if mutation:
            driver_mutation_annotation[gene] = mutation
        if indel:
            driver_indel_annotation[gene] = indel
        # TODO copy number
        # TODO fusion
    if v_type == "mutation":
        return driver_mutation_annotation
    elif v_type == "indel":
        return driver_indel_annotation
    else:
        print("Provided wrong variant type. Terminating ...")
        sys.exit()  # consider not terminating


def get_mechanistic_drugs(dataframe, knowledgebase):
    """query with hgnc id, select is_cancer_drug=true cases. Returns all the drugs mechanistically 
    targeting the affected gene."""
    gene_list = get_gene_list(dataframe)
    mechanistic_drug_targets = {}
    for gene in gene_list:
        mechanistic_drug = get_inner_dict(
            gene, "mechanistic_therapeutics_annotation", knowledgebase)
        if not mechanistic_drug:
            continue
        mechanistic_cancer_drug = [ # first filtering condiitions
            drug for drug in mechanistic_drug if (drug["is_cancer_drug"] == True and drug["interaction_type"] == "target")]
        if not mechanistic_cancer_drug:
            continue
        mechanistic_drug_targets[gene] = mechanistic_cancer_drug
    return mechanistic_drug_targets




def ann_pharm_variant(dataframe, genome, knowledgebase):
    """query with chr, start, stop, ref, alt, genome assembly version. Returns all the drugs 
    targeting the observd variant. """
    all_direct_targets = {}
    rows = get_variant_list(dataframe)
    for r in rows:
        direct_target_list = []
        gene = r[0]
        search_subset = {"chromosome": r[1], "start": r[2], "stop": r[3], "reference_base": r[4], "alteration_base": r[5], "assembly_version": genome}
        superset = get_inner_dict(gene, "variant_annotation", knowledgebase)
        if not superset:
            continue
        for variant in superset:
            v_index = superset.index(variant)
            search_set = superset[v_index]
            coverage = all(item in search_set.items() for item in search_subset.items())
            if coverage:
                direct_target_list.append(search_set)
        if not direct_target_list:
            continue
        if gene in all_direct_targets:
            all_direct_targets[gene].append(direct_target_list)
        else:
            all_direct_targets[gene] = direct_target_list
    return all_direct_targets




def ann_pharm_gene(dataframe, knowledgebase, genome):
    """query with hgnc id, returns all the drugs targeting the affected gene."""
    drug_targets_dict = {}
    gene_list = get_gene_list(dataframe)
    for gene in gene_list:
        drug_targets = get_inner_dict(gene, "variant_annotation", knowledgebase)
        if drug_targets:
            drug_targets_dict[gene] = drug_targets
    filtered = filter_genome_assembly(drug_targets_dict, genome)
    return filtered


def filter_variant_class(dictionary, variant_class_list):
    """Function to filter pharmaco-dyamics results based on variant class."""
    filtered = {}
    for key, value in dictionary.items():
        filtered[key] = [
            variant for variant in value if variant["variant_class"] in variant_class_list]
    return filtered

def filter_genome_assembly(dictionary, genome):
    """Function to filter pharmacogenomics dictionary(indirect one) based on the genome assembly"""
    filtered = {}
    for key, value in dictionary.items():
        filtered[key] = [
            variant for variant in value if variant["assembly_version"]==genome]
    return filtered
