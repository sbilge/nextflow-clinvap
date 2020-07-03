# Functions not specific to type (driver, mechanistic etc) coded here. 


# function to remove nested part of the dictionary given the key and the doctionary -
# aim is to remove the tumor_type ater creating the tumor string
# this can be used for more than one kind i.e for both drivers and mechanistic ones. That makes it a helper. move it to helpers.

#  input is either the list of keys or one key which is given as string? 
def remove_key_value(dictionary, *keys):
    """Delete key, value pairs from a dictionary."""
    for key in keys:
        if key in dictionary:
            del dictionary[key]
    return dictionary



def populate_list(dictionary, key, list_to_populate):
    """Functions to populate a list with valuesof the given key"""
    value = dictionary[key]
    if value != "null":
        list_to_populate.append(value)
    return list_to_populate


def populate_set(dictionary, key, set_to_populate):
    """Functions to populate a set with key values of the given key"""
    value = dictionary[key]
    set_to_populate.add(value)
    return set_to_populate


def create_tumor_string(tumor_dict_list):
    """Function to create tumor list string from tumor_type list of dicts"""
    if not tumor_dict_list:
        return None
    tumors = []
    for tumor in tumor_dict_list:
        abberration = tumor["abbr"]
        name = tumor["tumor_name"]
        if abberration != "":
            tumors.append(abberration)
        elif abberration == "" and name != "":
            tumors.append(name)
        else:
            continue
    tumor_string = "|".join(tumors)
    return tumor_string


def add_tumor_string(key, dictionary):
    """Function to add tumor abbreviation list to driver gene query resulting dict and
    removing tumor_type key, value from it."""
    annotation = dictionary[key]
    for info in annotation:
        info_index = annotation.index(info)
        tumor_info = info.get("tumor_type")
        tumor_string = create_tumor_string(tumor_info)
        if not tumor_string:
            tumor_string = ""
        if tumor_string == "|":
            tumor_string = ""
        annotation[info_index]["tumor_list"] = tumor_string
        remove_key_value(info, "tumor_type", "evidence_statement", "rating")
    return dictionary



   

