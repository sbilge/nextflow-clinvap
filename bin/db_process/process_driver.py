import db_process.helpers as helper
import pandas as pd
# For SNVs needs to be generalized. 

# FUNCTION FOR DRIVER GENES

# Driver gene dataframe with columns:
# GeneName
# Mutation conclusion?? 
# Observed mutation in that gene (missense, CNA, FUSION)
# Driver type
# Observed in cancer (abbreviations, DO name, DOID, ICD10 ID, only abbreviations will be used in the report template)
# Confidence - will I come uop with a better and proper confidence? or will I just remove this - how many sources has it been seen, is there a consensus for its role type, has it been observed in a cancer study
# References - thi will stay as it is. 



# MEDIUM TODO: SCRIPTTE BU YAPTIGIN MAPPING I ANLATMAN LAZIM

def assign_driver_role(role_set):
    """Function to collapse driver info from different sources into one, take driver role set, 
    assign driver role to a string"""
    if len(list(role_set)) == 1:
        driver_role = list(role_set)[0]
    elif role_set == {"TSG", "Oncogene"}:
        driver_role = "Oncogene/TSG"
    elif role_set == {"Oncogene/TSG", "TSG", "Oncogene"}:
        driver_role = "Oncogene/TSG"
    elif role_set == {"Oncogene/TSG", "Unknown", "TSG", "Oncogene"}:
        driver_role = "Oncogene/TSG"
    elif role_set == {"Unknown", "TSG", "Oncogene"}:
        driver_role = "Oncogene/TSG"
    elif role_set == {"Unknown", "TSG"}:
        driver_role = "TSG"
    elif role_set == {"Unknown", "Oncogene"}:
        driver_role = "Oncogene"
    elif role_set == {"Oncogene/TSG", 'Unknown'}:
        driver_role = "Oncogene/TSG"
    elif role_set == {"Oncogene/TSG", "TSG"}:
        driver_role = "Oncogene/TSG"
    elif role_set == {"Oncogene/TSG", "Unknown", "TSG"}:
        driver_role = "Oncogene/TSG"
    elif role_set == {"Oncogene/TSG", "Oncogene"}:
        driver_role == "Oncogene/TSG"
    return driver_role

    
def dict_to_dataframe(driver_dict):
    """Function to convert driver information dictionary into dataframe."""
    df_driver = pd.DataFrame.from_dict(
        driver_dict, orient="index").rename_axis("hgnc_id").reset_index()
    df_driver = df_driver.drop_duplicates().reset_index()
    return df_driver


