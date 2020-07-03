# VALID ONLY FOR SNVs

import table_content.handling as handle
from cyvcf2 import VCF

import sys
import pandas as pd



def read_vcf(inputfile):
    """Read vcf"""
    vcf = VCF(inputfile)
    return vcf


def get_fields(inputfile):
    """Get annotation field names - CSQ comes from Ensembl VEP. It might not work for SV or CNV"""
    vcf = read_vcf(inputfile)
    try:
        header = vcf.get_header_type(
            "CSQ")["Description"][51:].strip('"').split("|")
        return header
    except KeyError:
        print("Input file is not annotated. Run Ensembl VEP on input file. Terminating...")
        sys.exit()


def get_annotation(variant, header):
    """Get variant annotation info"""
    annotation = variant.INFO.get('CSQ').split('|')
    base = annotation[0]
    rest = annotation[1:]
    list_annotation = [rest[x: x + len(header) - 1]
                       for x in range(0, len(rest), len(header) - 1)]
    [ls.insert(0, base) for ls in list_annotation]
    return list_annotation


# Read VCF file and create dataframe minium variant level data (MVLD) according to
# Ritter et al. (https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-016-0367-z)

def parse_vcf(inputfile):
    """Parse and process VCF, create dataframe"""
    complete_line_list = list()
    fields = get_fields(inputfile)  # Get fields of VCF
    vcf = read_vcf(inputfile)
    for variant in vcf:
        chrom = variant.CHROM.replace("chr", "")
        if chrom != "X":
            chrom = int(chrom)
        start = variant.POS
        stop = variant.end
        vid = variant.ID
        ref = variant.REF
        qual = variant.QUAL
        vfilter = variant.FILTER
        alt = variant.ALT
        if len(alt) == 1:
            alt = alt[0]
        location = "{}:{}_{}/{}".format(chrom, start, ref, alt)
        line = [chrom, start, stop, location, vid, ref, alt, qual, vfilter]
        annotation = get_annotation(variant, fields)
        for ann in annotation:
            for ln in line:
                ann.append(ln)
            complete_line_list.append(ann)
    header = fields + ["chr", "start", "end", "location",
                       "vid", "ref", "alt", "qual", "filter"]  # create VCF dataframe
    dataframe = pd.DataFrame(complete_line_list, columns=header)

    dataframe["dbSNP"] = dataframe["Existing_variation"].str.findall("(rs.\d+)").apply(",".join)
    dataframe["COSMIC"] = dataframe["Existing_variation"].str.findall("(COSM.\d+)").apply(",".join)

    if dataframe.empty:
        handle.empty_mvld()
    else:
        return dataframe


def get_high_moderate_effect(dataframe):
    """Filter mvld, create mvld dataframe for high and moderate variant effects, works for Ensembl VEP output"""
    dataframe = dataframe[dataframe["PICK"] == "1"]

    if len(dataframe.index != 0):
        first_cond = dataframe[(dataframe["IMPACT"] == "HIGH") | (
            dataframe["IMPACT"] == "MODERATE")]
        second_cond = first_cond[~((first_cond["SIFT"] == "tolerated") & (
            first_cond["PolyPhen"] == "benign"))]
        third_cond = second_cond[~((second_cond["SIFT"] == "tolerated_low_confidence") & (
            second_cond["PolyPhen"] == "benign"))]
        return third_cond
    else:
        print("No variants found that passed default filters.")
        return pd.DataFrame()


def get_modifier_effect(dataframe):
    """Create mvld dataframe for variants that has modifying effect, works for Ensembl VEP output"""
    if len(dataframe.index != 0):
        cond = dataframe[(dataframe["IMPACT"] == "MODIFIER")]
        return cond
    else:
        print("No variants found that passed default filters.")
        return pd.DataFrame()


ONE_LETTER = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V", "Asx": "B", "Glx": "Z", "Xaa": "X"}


def one_letter_repr(dataframe):
    """Function to represent HGVSp with one letter aminoacids"""
    dataframe["aa"] = dataframe["HGVSp"].str.replace("ENSP.*p\.", "")
    dataframe["aa_1"] = dataframe["aa"].str.extract("([a-zA-Z]{3}(?=\d))")
    dataframe["aa_number"] = dataframe["aa"].str.extract( "((?<=[a-zA-Z]{3}).*(?=[a-zA-Z]{3}))")
    dataframe["aa_2"] = dataframe["aa"].str.extract("((?<=\d)[a-zA-Z]{3})")
    dataframe["aa_1"] = dataframe["aa_1"].map(ONE_LETTER)
    dataframe["aa_2"] = dataframe["aa_2"].map(ONE_LETTER)
    cols = ["aa_1", "aa_number", "aa_2"]
    dataframe["one_letter_repr"] = dataframe[cols].apply(lambda row: "".join(row.values.astype(str)), axis=1)
    dataframe = dataframe.drop(["aa", "aa_1", "aa_2", "aa_number"], axis = 1)
    return dataframe
    



