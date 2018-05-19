# Konner Kohl
# BI-GY 7663
# Final Project

# Import needed modules
import re
import os
from bioservices import KEGG
import csv

###### Database 1 ######

# KEGG

# Function to use KEGG to identify Drugs that target Colorectal Cancer
# Function will output drug names (keys) and IDs (values) as dictionary

def drug_dict(disease):

    # Create KEGG Object
    k = KEGG(verbose=False)

    # Create object for disease file
    dis = k.get(disease)

    # create dictionary of k.get() output with k.parse()
    # this is an extension of the KEGG class
    d = k.parse(dis)

    # Pull out Therapeutic drug information
    treatment_drugs = d["DRUG"]

    # Return dictionary of drugs
    return treatment_drugs


# Function to use drug IDs to gather target gene and pathway information
# outputs dictionary Gene names(keys) and HSA#(values)
# outputs dictionary Gene names(keys) and Drug IDs(values)


def drug_targets(drug_dic):
    # Create KEGG Object
    k = KEGG(verbose=False)

    # create empty list for drug IDs
    id_list = []

    # Create empty dictionary do add gene information to
    target_gene_dic = {}

    # create dictionary to link gene(key) and theraputic drug(value)
    gene_drug = {}

    # locate each drug id and add to list
    for value in drug_dic.values():
        id = re.findall(r"(D\d{5})", str(value))
        id_list.append(id[0])

    # Loop through drug IDs to gather information
    for drug_ID in id_list:

        # create object for drug information
        page = k.get(drug_ID)

        # create dictionary of drug information to isolate target information
        d = k.parse(page)

        # check for presence of target information
        if "TARGET" in d.keys():

            # isolate target information
            targ = d["TARGET"]

            # Remove pathways
            no_paths_pre = targ.split("  PATHWAY")

            # count spaces to identify presence of info
            spaces = targ.count(" ")

            # create list of genes
            gene_list = no_paths_pre[0].split("\n            ")

            # follow this if pathway section is present
            if spaces > 0:

                # loop through gene list
                for x in gene_list:
                    # separate gene names and HSA ID's
                    gene_split = x.split(" [")

                    # remove extras from gene name
                    y_split = gene_split[0].split(" ")

                    # add gene information to output dictionary
                    target_gene_dic[y_split[0]] = gene_split[1].strip("]")

                    # add gene and drug to output dictionary
                    gene_drug.setdefault(y_split[0], []).append(drug_ID)

            # if Gene doesn't have HSA# enter no value
            # also add gene to drug output dictionary
            else:
                target_gene_dic[no_paths_pre[0]] = ""

                for x in gene_list:
                    # separate gene names and HSA ID's
                    gene_split = x.split(" [")

                    # remove extras from gene name
                    y_split = gene_split[0].split(" ")

                    # add gene and drug to output dictionary
                    gene_drug.setdefault(y_split[0], []).append(drug_ID)

        else:
            pass

    return target_gene_dic, gene_drug


# Function to use gene HSA values to generate pathway information
# outputs dictionary Gene names(keys) and Pathways(values)

def target_paths(target_dict):
    # Create KEGG Object
    k = KEGG(verbose=False)

    # Create empty dictionary to output information
    gene_path = {}

    # start iterator
    i = 0

    # create list of targets
    target_names = list(target_dict.keys())

    # Loop through genes
    for HSA in target_dict.values():

        # Only use data where available
        if len(HSA) > 1:

            # get gene KEGG page
            page = k.get(HSA.lower())

            # isolate pathway information
            d = k.parse(page)

            # write pathway information to output dictionary
            if "PATHWAY" in d.keys():

                # create variable for pathways
                paths = d["PATHWAY"]

                # add pathway ids as list to gene name key
                gene_path[target_names[i]] = list(paths.keys())

                # increase iterator
                i += 1

            # add null value for no pathways
            else:
                gene_path[target_names[i]] = " "

                # increase iterator
                i += 1

        # Skip null values
        else:
            gene_path[target_names[i]] = " "

            # increase iterator
            i += 1

    return gene_path


###### End of Database 1 ######


###### Database 2 ######

# Cosmic

# Function will create dictionary of Gene names(keys) and drug resistance(values)
def drug_resistance(genes_dict, COSMIC_file):

    # create empty dictionary to write information to
    resistance_dict = {}

    # initialize dictionary with key values
    for y in genes_dict.keys():
        resistance_dict.setdefault(y, [None])

    # Open cosmic file
    with open(COSMIC_file, "r") as cosmic:

        # create object for cosmic file
        reader = csv.DictReader(cosmic, delimiter="\t")

        # loop through file by row
        for row in reader:

            # Only look at rows for genes of interest
            if row["Gene Name"] in genes_dict.keys() and row["Drug Name"] not in resistance_dict[row["Gene Name"]]:

                # create dictionary for drugs gene is resistant to
                resistance_dict.setdefault(row["Gene Name"], []).append(row["Drug Name"])

                # remove the None value from dictionary values after value is added
                if None in resistance_dict[row["Gene Name"]]:
                    resistance_dict[row["Gene Name"]].remove(None)

                else:
                    pass

            else:
                pass

    return resistance_dict

###### End Database 2 ######


###### Database 3 ######

# Function will parse GDC file for desired columns relating to targeted genes
# this function requires 6 inputs
# (1) The path to the GDC file, (2) pathway info, (3) resistance info, (4) gene ids, (5) thereputic drugs, (6) save location
def gdc_column_parse(file_path, paths_dic, resistance_dic, gene_dict, therep_drugs, save_path):

    # open GDC file
    with open(file_path, "r") as GDC:

        # Create output file
        gdc_parsed = open(save_path, "w")
        # create reader object
        reader = csv.DictReader(GDC, delimiter="\t")

        gdc_parsed.write("Hugo_Symbol" + "\t" + "KEGG_GENE" + "\t" + "Therapeutic_Drug" + "\t" "Drug Resistance" + "\t" +
                         "Pathways" + "\t" + "Entrez_Gene_Id" + "\t" + "Chromosome" + "\t" +
                         "Variant_Classification" + "\t" + "Variant_Type" + "\t" + "Tumor_Seq_Allele1" + "\t" +
                         "Tumor_Seq_Allele2" + "\t" + "dbSNP_RS" + "\t" + "Tumor_Validation_Allele1" + "\t" +
                         "Tumor_Validation_Allele2" + "\t" + "all_effects" + "\t" + "Allele" + "\t" +
                         "Gene" + "\t" + "RefSeq" + "\t" + "SIFT" + "\t" + "PolyPhen" + "\t" + "IMPACT" + "VARIANT_CLASS" + "\n")

        # only use columns of desired genes
        for row in reader:

            if row["Hugo_Symbol"] in paths_dic.keys():
                gdc_parsed.write(row["Hugo_Symbol"] + "\t" + str(gene_dict[row["Hugo_Symbol"]]) + "\t" +
                                 str(therep_drugs[row["Hugo_Symbol"]]) + "\t" + str(resistance_dic[row["Hugo_Symbol"]]) + "\t" +
                                 str(paths_dic[row["Hugo_Symbol"]]) + "\t" + row["Entrez_Gene_Id"] + "\t" +
                                 row["Chromosome"] + "\t" + row["Variant_Classification"] + "\t" +
                                 row["Variant_Type"] + "\t" + row["Tumor_Seq_Allele1"] + "\t" +
                                 row["Tumor_Seq_Allele2"] + "\t" + row["dbSNP_RS"] + "\t" +
                                 row["Tumor_Validation_Allele1"] + "\t" + row["Tumor_Validation_Allele2"] + "\t" +
                                 row["all_effects"] + "\t" + row["Allele"] + "\t" + row["Gene"] + "\t" +
                                 row["RefSeq"] + "\t" + row["SIFT"] + "\t" + row["PolyPhen"] + "\t" + row["IMPACT"] +
                                 row["VARIANT_CLASS"] + "\n")
        gdc_parsed.close()

    return

###### End of Functions ######


# Create variable for Disease ID
colorectal_cancer = "H00020"

# Enter GDC file location
GDC_file = "/Users/konnerkohl/BI-GY_7663_Final/TCGA.COAD.mutect.txt"

# Enter COSMIC file Location
COSMIC_file = "/Users/konnerkohl/BI-GY_7663_Final/CosmicResistanceMutations.tsv"

# Create variable for current working directory
cwd = os.getcwd()

# create save path for output file
file_path = cwd + "/" + str(colorectal_cancer) + "_Annotated.txt"

# Run through each function
func1_out = drug_dict(colorectal_cancer)
func2_out = drug_targets(func1_out)
func3_out = target_paths(func2_out[0])
func4_out = drug_resistance(func2_out[0], COSMIC_file)
gdc_column_parse(GDC_file, func3_out, func4_out, func2_out[0], func2_out[1], file_path)

