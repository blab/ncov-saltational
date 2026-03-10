import jsonlines
import argparse
import sys
import pandas as pd
import gzip

if __name__ == "__main__":

    # Load arguments
    parser = argparse.ArgumentParser(
        description="Get metadata file (stored in the Taxonium JSONL) and output as TSV",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--input_taxonium', type=str, required=True, help='Taxonium JSONL file')
    parser.add_argument('--output_metadata', type=str, required=True, help='TSV metadata')
    parser.add_argument('--name_reference', type=str, required=True, help='Reference sequence name')
    parser.add_argument('--date_reference', type=str, required=True, help='Reference sequence date')
    parser.add_argument('--country_reference', type=str, required=True, help='Reference sequence country')
    parser.add_argument('--pango_reference', type=str, required=True, help='Reference sequence pango lineage')
    args = parser.parse_args()

    
    # Initate metadata with reference sequence information
    rows_viridian = []
    rows_viridian.append([args.name_reference, args.country_reference, args.date_reference, args.pango_reference, ""])

    # Read Taxonium tree
    with gzip.open(args.input_taxonium, 'rt') as f:
        reader = jsonlines.Reader(f)
        for obj in reader:
            if obj.get('name') != None:
                rows_viridian.append([obj.get('name'), obj.get('meta_Country'), obj.get('meta_Date_tree'), obj.get('meta_Pangolin_usher'),  obj.get('meta_Platform')])

    # Convert to dataframe
    df_viridian = pd.DataFrame(rows_viridian)

    # Name columns
    df_viridian = df_viridian.rename(columns={
        0: 'strain',
        1: 'country',
        2: 'date',
        3: 'pangolin_usher',
        4: 'platform'
    })

    # Remove internal nodes
    metadata_samples = df_viridian[~df_viridian['strain'].str.startswith('node_')]


    
    # Write results to output file
    metadata_samples.to_csv(args.output_metadata, sep = "\t", index = False)
