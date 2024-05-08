import pandas as pd
import numpy as np
import requests
#from pyspark.sql import SparkSession
from extract import get_raw_data, get_protein_location, get_protein_function, get_protein_process, get_mapped_ids, get_interactions

#1. Transform raw data from proteinGroups.txt file, extract protein ids and gene names for use. Create dataframe

def tf_raw_data():
    raw_data_df = get_raw_data()
    raw_data_df['Gene names'] = raw_data_df['Gene names'].astype(str) #Bug in the source data - wrong data type for gene name (float, not string)
    protein_gene_df0 = raw_data_df[['Protein IDs','Gene names']]

    #Quality check data, remove chaff and duplicates
    protein_gene_df0.dropna(subset=['Gene names'], inplace=True) #drop NaN values from df
    protein_gene_df0['Gene names'] = [i.upper().strip().split(';')[0] for i in protein_gene_df0['Gene names']]
    protein_gene_df0['Protein IDs'] = [i.upper().strip().split(';')[0] for i in protein_gene_df0['Protein IDs']]
    protein_gene_df1 = protein_gene_df0.drop_duplicates(subset='Protein IDs') #Only keep rows which have unique values in 'Protein IDs' column
    protein_gene_df2 = protein_gene_df1[~protein_gene_df1['Protein IDs'].str.startswith(('REV', 'CON'))]     #Remove outliers which remain: thos ewhich start with 'REV" or 'CON' (contaminant)
    protein_gene_df2['Gene names'].replace('nan', np.nan, inplace=True) #replace 'nan' strings with NaN type
    protein_gene_df2.dropna(subset=['Gene names'], inplace=True) #drops any na values in both columns
    protein_gene_dict = dict(zip(protein_gene_df2.iloc[47:111, 0], protein_gene_df2.iloc[47:111, 1])) # {protein_id : gene name} - examples for demo
    protein_gene_df_final = protein_gene_df2
    protein_gene_df_final.reset_index(drop=True, inplace=True) #reset index post-transformation

    return protein_gene_dict, protein_gene_df_final


#2. Extract protein location information and create dataframe

def tf_protein_locations():
    protein_gene_dict, _ = tf_raw_data()
    protein_ids = protein_gene_dict.keys()
    list_of_locations_dict = [{'Protein IDs': protein_id, 'Location': ', '.join(get_protein_location(protein_id).get(protein_id, []))} for protein_id in protein_ids]
    list_of_locations_dict = [value for value in list_of_locations_dict if value['Location']] #removes any dicts from list with no value pairing
    locations_df = pd.DataFrame(list_of_locations_dict)
    return locations_df


#3. Extract protein function information and create dataframe

def tf_protein_functions():
    protein_gene_dict, _ = tf_raw_data()
    protein_ids = protein_gene_dict.keys()
    list_of_functions_dict = [get_protein_function(protein_id) for protein_id in protein_ids]    
    unpack_functions_dict = [{'Protein IDs': protein_id, 'Function': function} for element in list_of_functions_dict for protein_id, function in element.items()]
    unpack_functions_dict = [value for value in unpack_functions_dict if value['Function']]
    functions_df = pd.DataFrame(unpack_functions_dict)
    return functions_df


#4. Extract protein process information and create dataframe

def tf_protein_processes():
    protein_gene_dict, _ = tf_raw_data()
    protein_ids = protein_gene_dict.keys()
    list_of_process_dict = [get_protein_process(protein_id) for protein_id in protein_ids]    
    unpack_process_dict = [{'Protein IDs': protein_id, 'Process': process} for element in list_of_process_dict for protein_id, process in element.items()]
    unpack_process_dict = [value for value in unpack_process_dict if value['Process']]
    process_df = pd.DataFrame(unpack_process_dict)
    return process_df


#5. Use UniProt to STRING protein ID mappings to extract interaction data from STRING database. Inner join based gene name (creates merged_df)

def tf_interactions():
    protein_gene_dict, protein_gene_df_final = tf_raw_data()
    protein_ids = protein_gene_dict.keys()
    gene_names = protein_gene_dict.values()
    Uniprot_to_STRING_mappings = get_mapped_ids(gene_names) #obtain id-to-id mappings (UniProt - String)
    STRING_interaction_dicts = get_interactions(gene_names, protein_ids)
    interactions_df = pd.DataFrame(STRING_interaction_dicts)
    interactions_df = interactions_df.drop_duplicates(subset=['query_protein_STRING_id', 'interacting_protein_STRING_id'], ignore_index=True)
    merged_df = protein_gene_df_final.merge(interactions_df, how='inner', left_on='Gene names', right_on = 'query_gene_name') #inner join based on shared gene_name values
    pd.set_option('display.max_columns', 40)
    return interactions_df, merged_df


# TEST ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#1.

# protein_gene_dict, protein_gene_df_final = tf_raw_data()
# print(protein_gene_df_final)

#2.

# locations_df = tf_protein_locations()
# print(locations_df)

#3. 

# pd.set_option('display.max_rows', 400)
# functions_df = tf_protein_functions()
# print(functions_df)

#4. 

# process_df = tf_protein_processes()
# print(process_df)

#5.

# interactions_df, merged_df= tf_interactions()#, merged_df = 
# print(merged_df)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



if __name__ == "__main__":
    pass