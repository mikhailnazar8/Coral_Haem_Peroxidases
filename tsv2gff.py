import pandas as pd
import argparse

def convert_to_gff(input_file, output_file):
    # Read the input file into a pandas DataFrame
    df = pd.read_csv(input_file, sep='\t')
    
    # Define GFF columns and prepare for GFF format
    gff_columns = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
    gff_data = []
    
    for _, row in df.iterrows():
        # Extract core GFF columns
        seqid = row['id']
        source = "interproscan"  
        type_ = "CDS"    
        start = row['start']
        end = row['end']
        score = row['score'] if pd.notna(row['score']) else '.'
        strand = "."  # Strand information is missing in data, so using '.'
        phase = "."    # Phase is not applicable in this context, so using '.'
        
        # Extract and format attributes
        attributes = (
            f"ID={row['id']};"
            f"MD5={row['md5']};"
            f"Length={row['length']};"
            f"Analysis={row['analysis']};"
            f"Sig_Acc={row['sig_acc']};"
            f"Sig_Desc={row['sig_desc'].replace(' ', '_')};"
            f"Status={row['status']};"
            f"IPR_Acc={row['ipr_acc']};"
            f"IPR_Desc={row['ipr_desc'].replace(' ', '_')};"
            f"GO_Term={row['goterm']}"
        )
        
        # Append the row to gff_data
        gff_data.append([seqid, source, type_, start, end, score, strand, phase, attributes])
    
    # Create a DataFrame for GFF data and save it as a GFF file
    gff_df = pd.DataFrame(gff_data, columns=gff_columns)
    gff_df.to_csv(output_file, sep='\t', index=False, header=False, quoting=3)

def main():
    # Argument parser setup
    parser = argparse.ArgumentParser(description='Convert TSV formatted data to GFF format.')
    parser.add_argument('input_file', metavar='<TSV FILE>', type=str, help='Input TSV file')
    parser.add_argument('output_file', metavar='<GFF FILE>', type=str, help='Output GFF file')

    # Parse arguments
    args = parser.parse_args()

    # Convert the TSV file to GFF
    convert_to_gff(args.input_file, args.output_file)

if __name__ == "__main__":
    main()
