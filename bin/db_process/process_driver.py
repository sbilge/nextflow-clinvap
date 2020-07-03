import db_process.helpers as helper
import pandas as pd


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


