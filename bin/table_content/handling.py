import sys
import pandas as pd


def empty_mvld():
    """Function to exit program if the main ataframe calculated from annotated VCF empty."""
    message = "Input VFC is empty. Terminating..."
    print(message)
    sys.exit()


def empty_database_result():
    """Function to exit program if there are not database results returned for the input VCF."""
    message = "No database result. Terminating..."
    print(message)
    sys.exit()

def empty_high_impact_mvld():
    """Function to print the wrning message of the hogha nad moderate effect mvld table is empty."""
    message = "MVLD does not contain variants with high or moderate impact. No variants found that passed default filters. Modifier effect variants will be used for further analyis, if it is not empty."
    print(message)


def empty_driver_annotation():
    """Function to print a warning if querying the knowledgebase for driver genes returns empty. Sets driver flag to false"""
    message = "No driver gene identified. Driver gene table will be empty."
    print(message)


def empty_mechanistic():
    """Function to print warning that there is no hit for mechanistic drug targets."""
    message = "No mechanistic drugs targeting distrupted genes was identified. Mechanistic drug targets table will be empty."
    print(message)


def empty_approved_mechanistic():
    """Function to print warning when there are not approved drugs in mechanistic drug list"""
    message = "No approved drugs was returned for the distrupted genes. Investigational drugs will be used for further analysis, if it is not empty."
    print(message)


def empty_direct_variant_annotation():
    """Function to print a warning message if query to get variant annotation for observed variants is empty."""
    message = "No direct pharmacodynamics information was found for observed variants."
    print(message)


def empty_variant_annotation():
    """Function to print a warning message if query to get variant annotation for observed genes is empty."""
    message = "No pharmacodynamics information was found for observed genes."
    print(message)



def pharm_columns(combined=None):
    """Function to return empty dataframes with pharmacogenomics dataframe column names"""
    if combined:
        df = pd.DataFrame(columns=['combined_variant', 'drug_approval_status', 'drug_class',
                                   'drug_drug_relationship', 'drug_name', 'drugbank_id', 'evidence_level',
                                   'paired_gene', 'paired_gene_hgnc', 'paired_variant', 'reference_id',
                                   'reference_source', 'source_name', 'source_pmid',
                                   'variant_drug_association', 'tumor_list', 'hgnc_id', 'variant',
                                   'variant_class', 'variant_type', 'chromosome', 'assembly_version',
                                   'alteration_base', 'reference_base', 'start', 'stop'])
    else:
        df = pd.DataFrame(columns=['drug_approval_status', 'drug_class', 'drug_drug_relationship',
                                'drug_name', 'drugbank_id', 'evidence_level', 'reference_id',
                                'reference_source', 'source_name', 'source_pmid', 'targeting',
                                'variant_drug_association', 'tumor_list', 'hgnc_id', 'variant',
                                'variant_class', 'variant_type', 'chromosome', 'assembly_version',
                                'alteration_base', 'reference_base', 'start', 'stop'])
    return df

def empty_dataframe_direct_pharmacogenomics():
    """Function to return empty dataframe with column names if direct pharmacogenomics dataframe is empty."""
    message = "No pharmacogenomics information found for observed variants."
    print(message)
    place_holder = pharm_columns()
    return place_holder


def empty_dataframe_direct_pharmacogenomics_combined_variants():
    """Not implemented. Always empty for MUT, MUT"""
    pass



def empty_dataframe_pharmacogenomics():
    """Function to return empty dataframe with column names if pharmacogenomics dataframe is empty."""
    message = "No pharmacogenomics information found for observed genes."
    print(message)
    place_holder = pharm_columns()
    return place_holder


def empty_dataframe_pharmacogenomics_combined():
    """Function to return empty dataframe with column names if pharmacogenomics combined variants dataframe 
    is empty."""
    message = "No pharmacogenomics information found for observed gene combinations."
    print(message)
    place_holder = pharm_columns(combined = 1)
    return place_holder


def adverse_effect_columns(combined=None):
    if combined:
        df = pd.DataFrame(columns=['combined_variant', 'drug_approval_status', 'drug_class', 'drug_name',
                                   'drugbank_id', 'evidence_level', 'paired_gene', 'paired_gene_hgnc',
                                   'paired_variant', 'reference_id', 'reference_source', 'source_name',
                                   'source_pmid', 'tumor_type', 'variant_drug_association', 'hgnc_id',
                                   'variant', 'variant_class', 'variant_type', 'chromosome',
                                   'assembly_version', 'alteration_base', 'reference_base', 'start',
                                   'stop'])
    else:
        df = pd.DataFrame(columns=['drug_approval_status', 'drug_class', 'drug_name', 'drugbank_id',
                                   'evidence_level', 'evidence_statement', 'reference_id',
                                   'reference_source', 'source_name', 'source_pmid', 'tumor_type',
                                   'variant_drug_association', 'hgnc_id', 'variant', 'variant_class',
                                   'variant_type', 'chromosome', 'assembly_version', 'alteration_base',
                                   'reference_base', 'start', 'stop'])
    return df


def empty_dataframe_direct_adverse():
    """Function to return empty dataframe with column names if direct adverse effect dataframe is empty."""
    message = "No adverse efect information was found for observed variants."
    print(message)
    place_holder = adverse_effect_columns()
    return place_holder


def empty_dataframe_adverse():
    """Function to return empty dataframe with column names if adverse effect dataframe is empty."""
    message = "No adverse effect information was found for observed genes."
    print(message)
    place_holder = adverse_effect_columns()
    return place_holder


def empty_dataframe_adverse_combined():
    """Function to return empty dataframe with column names if adverse effect combined variants 
    dataframe is empty."""
    message = "No adverse effect information was found for observed gene combinations."
    print(message)
    place_holder = adverse_effect_columns(combined=1)
    return place_holder



def empty_reference_place_holder():
    """Function to return empty reference dataframe with column names if variant annotation 
    (both for genes and variants) returns empty."""
    df = pd.DataFrame(columns=['reference_id', 'reference_source'])
    return df


def empty_pharmacodynamics_place_holder(adverse=None):
    """Function to return empty pharmacodynamics dataframes with columns if variant annotation 
    (both for genes and variants) empty"""
    message = "Variant annotation is empty. Place holders will be returned."
    print(message)
    df = pd.DataFrame(columns=['SYMBOL', 'drug_name', 'variant_drug_association', 'tumor_list',
                               'hgnc_id', 'variant', 'variant_type', 'match_level', 'reference_id'])
    if adverse:
        df = df.drop(columns=['tumor_list'])
    return df

