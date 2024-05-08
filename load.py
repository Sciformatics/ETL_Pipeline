from sqlalchemy import create_engine #import sqlalchemy for database
from transform import tf_raw_data, tf_protein_locations, tf_protein_functions, tf_protein_processes, tf_interactions
import os
 

#instantiate the db engine in your pwd
pwd = os.getcwd()
engine = create_engine(f'sqlite:///{pwd}/test_database.db')

#Load transformed dataframes to the SQL db

def load_raw_data():
    _, protein_gene_df_final = tf_raw_data()
    protein_gene_df_final.to_sql('protein_gene_ids', engine, if_exists='replace')

def load_locations():
    locations_df = tf_protein_locations()
    locations_df.to_sql('protein_locations', engine, if_exists='replace')

def load_functions():
    functions_df = tf_protein_functions()
    functions_df.to_sql('protein_functions', engine, if_exists='replace')

def load_processes():
    process_df = tf_protein_processes()
    process_df.to_sql('protein_processes', engine, if_exists='replace')

def load_interactions():
    interactions_df, merged_df = tf_interactions()
    interactions_df.to_sql('protein_interactions', engine, if_exists='replace')
    merged_df.to_sql('merged_protein_interactions', engine, if_exists='replace')

if __name__ == "__main__":
    pass