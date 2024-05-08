#modules
import requests, os, warnings
import pandas as pd
import xml.etree.ElementTree as ET


#~~~~~~~~~~~~~~~~~~~ GET functions ~~~~~~~~~~~~~~~~~~~~~
warnings.filterwarnings("ignore") 

#1. Read, parse and extract data from proteinGroups file

def get_raw_data():
    pwd = os.getcwd()
    raw_data_df = pd.read_table(f'{pwd}/proteinGroups.txt', sep='\t') #read in raw data file
    return raw_data_df
 

#2. Retrieve location of protein in the cell

def get_protein_location(protein_id):
    url = f"https://www.uniprot.org/uniprot/{protein_id}.xml"

    response = requests.get(url)
    if response.ok:
        xml_data = response.text
        root = ET.fromstring(xml_data)
        ns = {"uniprot": "http://uniprot.org/uniprot"}

        #Finds all subcellular location annotations
        locations = root.findall(".//uniprot:subcellularLocation", namespaces=ns)

        subcellular_localisations = []
        for location in locations:
            localisation = location.find("uniprot:location", namespaces=ns).text #Extract the subcellular localisation
            subcellular_localisations.append(localisation)
        localisation_dict = {protein_id: subcellular_localisations}
        return localisation_dict
    

# 3. Retrieve protein function within the cell

def get_protein_function(protein_id):
    url = f"https://www.uniprot.org/uniprot/{protein_id}.xml"
    response = requests.get(url)
    function_dict= {}
    if response.ok:
        xml_data = response.text
        root = ET.fromstring(xml_data)
        ns = {"uniprot": "http://uniprot.org/uniprot"}

        #Find the protein function annotation
        function = root.find(".//uniprot:comment[@type='function']", namespaces=ns)
        if function is not None:
            #Extract the protein function descrption
            function_description = function.find("uniprot:text", namespaces=ns).text
            function_description = function_description.strip().split('.')[0] #Splits at '.'
            function_dict[protein_id] = function_description
    return function_dict
        

# 4. Retrieve information on where the protein operates/which process it is part of

def get_protein_process(protein_id):
    url = f"https://www.uniprot.org/uniprot/{protein_id}.xml"
    response = requests.get(url)
    process_dict = {}
    if response.ok:
        root = ET.fromstring(response.text)
        ns = {"uni": "http://uniprot.org/uniprot"}
        process_elements = root.findall(".//uni:dbReference[@type='GO']", namespaces=ns)
        for process_element in process_elements:
            process_id = process_element.get("id")
            process_term = process_element.find("./uni:property[@type='term']", namespaces=ns).get("value")
            process_term = process_term[2:] #remove 'P:' prefix to process description
            process_dict[protein_id] = process_term


    return process_dict

# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# #Note 1: We've now got the 100 unique protein ids (behaving as primary keys)
# #Note 2: We've obtained each protein's function, its location within the cell and the process it's involved in
# #Note 3: Now, we can see which of these 100 proteins interact with one another - time to map the protein_ids to the STRING database to \
# # retrieve interaction info :-)
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 5. Retrieve STRING protein ids for the same proteins, map to UniProt protein ids

#WORKS

def get_mapped_ids(gene_names):
    # protein_gene_dict = get_protein_ids()
    # protein_ids = protein_gene_dict.values()
    Uniprot_to_STRING_mappings = {}
    string_api_url = "https://version-11-5.string-db.org/api"
    output_format = "tsv-no-header"
    method = "get_string_ids"
    params = {
        "identifiers" : "\r".join(gene_names),
        "species" : 9606, #species NCBI identifier 
        "limit" : 1, #only one (best) identifier per input protein
        "echo_query" : 1, #see the input identifiers in the output
        "columns": "id, proteinNames",
    }

    request_url = "/".join([string_api_url, output_format, method])
    results = requests.post(request_url, data=params)

    for line in results.text.strip().split("\n"):
        l = line.split("\t")
        if len(l) >= 3:  #Ensures that the list has at least 3 elements
            uniprot_identifier = l[0]
            string_identifier = l[2].strip().split(';')[0]
            Uniprot_to_STRING_mappings[uniprot_identifier] = string_identifier
    return Uniprot_to_STRING_mappings


# 6. Retrieve the interaction information showing which of the proteins interact, including the db and combined score

def get_interactions(gene_names, protein_ids):
    STRING_interaction_dicts = []
    Uniprot_to_STRING_mappings = get_mapped_ids(protein_ids)
    #API call
    string_api_url = "https://version-11-5.string-db.org/api"
    output_format = "tsv-no-header"
    method = "network"
    request_url = "/".join([string_api_url, output_format, method])
    params = {
        "identifiers": "%0d".join(gene_names), # Use values of the mappings
        "species": 9606,
    }
    response = requests.post(request_url, data=params)
    for line in response.text.strip().split("\n"):
        interaction_data = line.strip().split("\t")
        interaction_dict = {
            'query_protein_STRING_id': interaction_data[0], #query protein id in STRING identifier format
            'query_gene_name': interaction_data[2], #query gene name
            'interacting_protein_STRING_id': interaction_data[1], #interacting protein id in STRING identifier format
            'interacting_gene_name': interaction_data[3], #interacting gene name
            'db_score': '{:.2f}'.format(float(interaction_data[11])), #database score
            'combined_score':'{:.2f}'.format(float(interaction_data[5])) #combined score
        }
        STRING_interaction_dicts.append(interaction_dict) #append as a list of dictionaries

    return STRING_interaction_dicts


if __name__ == "__main__":
    pass