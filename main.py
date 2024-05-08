from load import load_raw_data, load_locations, load_functions, load_processes, load_interactions, engine

def main():

    #1. Load transformed raw data table containing protein ids and gene names
    load_raw_data()

    #2. Load the protein location, function and process annotation data
    load_locations()
    load_functions()
    load_processes()

    #3. load the interaction data obtained from STRING protein interaction database
    load_interactions()

if __name__ == "__main__":
    main()
