import pandas as pd
import argparse

def parse_attributes(attributes_str):
    """
    Parse the GFF attributes column into a dictionary.
    """
    attributes = {}
    for attribute in attributes_str.split(';'):
        key_value = attribute.split('=')
        if len(key_value) == 2:
            key, value = key_value
            attributes[key] = value
    return attributes

def convert_to_tsv(gff_file, tsv_file):
    # Define the columns for the TSV file
    tsv_columns = ["id", "start", "end", "score", "md5", "length", "analysis", 
                   "sig_acc", "sig_desc", "status", "ipr_acc", "ipr_desc", "goterm"]
    
    tsv_data = []
    
    # Read the GFF file
    with open(gff_file, 'r') as gff:
        for line in gff:
            if line.startswith('#'):
                continue  # Skip comment lines
            
            parts = line.strip().split('\t')
            seqid = parts[0]
            start = int(parts[3])
            end = int(parts[4])
            score = parts[5] if parts[5] != '.' else None
            attributes_str = parts[8]
            
            # Parse the attributes column
            attributes = parse_attributes(attributes_str)
            
            # Extract attributes
            md5 = attributes.get("MD5", "")
            length = attributes.get("Length", "")
            analysis = attributes.get("Analysis", "")
            sig_acc = attributes.get("Sig_Acc", "")
            sig_desc = attributes.get("Sig_Desc", "").replace('_', ' ')
            status = attributes.get("Status", "")
            ipr_acc = attributes.get("IPR_Acc", "")
            ipr_desc = attributes.get("IPR_Desc", "").replace('_', ' ')
            goterm = attributes.get("GO_Term", "")
            
            # Append the row to tsv_data
            tsv_data.append([seqid, start, end, score, md5, length, analysis, 
                             sig_acc, sig_desc, status, ipr_acc, ipr_desc, goterm])
    
    # Create a DataFrame for TSV data and save it as a TSV file
    tsv_df = pd.DataFrame(tsv_data, columns=tsv_columns)
    tsv_df.to_csv(tsv_file, sep='\t', index=False)

def main():
    # Argument parser setup
    parser = argparse.ArgumentParser(description='Convert GFF formatted data back to TSV format.')
    parser.add_argument('gff_file', metavar='<GFF FILE>', type=str, help='Input GFF file')
    parser.add_argument('tsv_file', metavar='<TSV FILE>', type=str, help='Output TSV file')

    # Parse arguments
    args = parser.parse_args()

    # Convert the GFF file to TSV
    convert_to_tsv(args.gff_file, args.tsv_file)

if __name__ == "__main__":
    main()