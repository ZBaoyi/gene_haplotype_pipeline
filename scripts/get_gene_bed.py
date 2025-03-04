import argparse
import subprocess
import sys
import os

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Extract genomic regions from GFF based on gene ID and region type, then use PLINK to extract genotype data."
    )
    parser.add_argument('--gene', type=str, required=True,
                        help='Target gene ID, e.g., GeneX')
    parser.add_argument('--region', type=str, required=True, choices=['promoter', 'gene_body', 'promoter_gene_body'],
                        help='Region type: promoter, gene_body, or promoter_gene_body')
    parser.add_argument('--gff', type=str, required=True,
                        help='Path to the input GFF file')
    parser.add_argument('--plink_input', type=str, required=True,
                        help='PLINK input file (prefix or specific file, depending on input type)')
    parser.add_argument('--data_type', type=str, required=True, choices=['bfile', 'file', 'vcf'],
                        help='PLINK data type: bfile: bed/bim/fam (default), file: ped/map or vcf')
    parser.add_argument('--plink_out', type=str, required=True,
                        help='PLINK output prefix')
    return parser.parse_args()

def extract_gene_region(gff_file, gene_id, region_type):
    gene_info = None
    gene_start = None
    gene_end = None
    # Extract gene record from GFF (feature type "gene")
    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            if fields[2] == "mRNA" and f"ID={gene_id}" in fields[8]:
                gene_info = fields
                try:
                    gene_start = int(fields[3])
                    gene_end = int(fields[4])
                except ValueError:
                    sys.exit("Invalid coordinate format in gene record of GFF file.")
                break
    if gene_info is None:
        sys.exit(f"Gene ID {gene_id} not found in {gff_file}.")
    chrom = gene_info[0]
    if chrom.startswith("chr"):
        chrom = chrom[3:]
    strand = gene_info[6]
    # For plus strand, TSS is gene_start; for minus strand, TSS is gene_end.
    tss = gene_start if strand == '+' else gene_end

    # Determine region boundaries based on region type
    if region_type == 'promoter':
        if strand == '+':
            region_start = max(gene_start - 1000, 1)
            region_end = gene_start - 1
        else:
            region_start = gene_end + 1
            region_end = gene_end + 1000
    elif region_type == 'gene_body':
        region_start = gene_start
        region_end = gene_end
    elif region_type == 'promoter_gene_body':
        if strand == '+':
            region_start = max(gene_start - 1000, 1)
            region_end = gene_end
        else:
            region_start = gene_start
            region_end = gene_end + 1000

    # Extract all feature information for the gene from GFF (excluding the gene record itself)
    feature_list = []
    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            # Check if the record belongs to the target gene and is not a gene record
            if f"Parent={gene_id}" in fields[8] and fields[2] != "mRNA":
                
                try:
                    feat_start = int(fields[3])
                    feat_end = int(fields[4])
                except ValueError:
                    continue
                feat_type = fields[2]
                feature_list.append((feat_start, feat_end, feat_type))
    feature_list.sort(key=lambda x: x[0])
    return chrom, region_start, region_end, tss, strand, gene_start, gene_end, feature_list


def build_plink_command(data_type, plink_input, chrom, region_start, region_end, plink_out):
    # Select PLINK input type
    if data_type == 'bfile':
        base_cmd = ['plink', '--bfile', plink_input]
    elif data_type == 'ped':
        base_cmd = ['plink', '--file', plink_input]
    elif data_type == 'vcf':
        base_cmd = ['plink', '--vcf', plink_input]
    else:
        sys.exit("Unsupported data type. Choose from: bfile, ped, or vcf.")

    # Construct the PLINK command
    base_cmd.extend([
        '--chr', str(chrom),
        '--from-bp', str(region_start),
        '--to-bp', str(region_end),
        '--recode','--double-id', 
        '--out', plink_out
    ])
    return base_cmd

def output_gene_info_file(plink_out, tss, strand, gene_start, gene_end, feature_list):
    vcf_file = plink_out + ".map"
    info_file = plink_out + "_gene_info.txt"
    if not os.path.exists(vcf_file):
        sys.exit(f"map file {vcf_file} not found.")
    try:
        with open(vcf_file, 'r') as vf, open(info_file, 'w') as of:
            # Write header: CHROM, ABS_POS, REL_POS, REGION_CATEGORY
            of.write("CHROM\tABS_POS\tREL_POS\tREGION_CATEGORY\n")
            for line in vf:
                if line.startswith('#'):
                    continue
                fields = line.strip().split()
                chrom_val = fields[0]
                try:
                    abs_pos = int(fields[3])
                except ValueError:
                    continue
                # Compute relative position: for plus strand: abs_pos - tss + 1; for minus strand: tss - abs_pos + 1
                if strand == '+':
                    rel_pos = abs_pos - tss + 1
                else:
                    rel_pos = tss - abs_pos + 1
                # Determine region category
                # If relative position is negative, it is promoter (manual designation)
                if rel_pos < 0:
                    region_category = "promoter"
                else:
                    # For positions in gene body, check all features from GFF that cover the SNP
                    matched_features = []
                    for feat_start, feat_end, feat_type in feature_list:
                        if feat_start <= abs_pos <= feat_end:
                            matched_features.append(feat_type)
                    if matched_features:
                        region_category = ";".join(sorted(set(matched_features)))
                    else:
                        region_category = "intron"
                of.write(f"{chrom_val}\t{abs_pos}\t{rel_pos}\t{region_category}\n")
        print(f"Gene info mapping file written to {info_file}")
    except Exception as e:
        sys.exit(f"Error writing gene info file: {e}")

def main():
    args = parse_arguments()
    # Extract gene region and feature information from GFF
    chrom, region_start, region_end, tss, strand, gene_start, gene_end, feature_list = extract_gene_region(args.gff, args.gene, args.region)
    print(f"Extracted region for gene {args.gene}: {chrom}\t{region_start}\t{region_end}")
    plink_cmd = build_plink_command(args.data_type, args.plink_input, chrom, region_start, region_end, args.plink_out)
    print("Executing PLINK command:")
    print(" ".join(plink_cmd))
    try:
        result = subprocess.run(plink_cmd, check=True, capture_output=True, text=True)
        print("PLINK executed successfully.")
        print(result.stdout)
    except subprocess.CalledProcessError as e:
        sys.exit(f"PLINK execution failed: {e.stderr}")
    # Output gene info mapping file including absolute and relative positions and region category
    output_gene_info_file(args.plink_out, tss, strand, gene_start, gene_end, feature_list)

if __name__ == '__main__':
    main()