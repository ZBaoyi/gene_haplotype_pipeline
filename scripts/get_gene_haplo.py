import argparse
iupac_to_nucleotides = {
    'AG': "R",
    'GA': "R",
    'CT': "Y",
    'TC': "Y",
    'GC': "S",
    'CG': "S",
    'AT': "W",
    'TA': "W",
    'GT': "K",
    'TG': "K",
    'AC': "M",
    'CA': "M",
}

def parse_arguments():
    parser = argparse.ArgumentParser(description="Generate haplotype information from PLINK map/bed and gene info file.")
    parser.add_argument('--ped', type=str, required=True, help='Input PLINK .ped file')
    parser.add_argument('--output', type=str, default="haplotype_summary.txt", help='Output file for haplotype summary')
    return parser.parse_args()


def generate_haplotypes(ped_file):
    hap_dic = {}
    with open(ped_file, 'r') as f:
        for line in f:
            l = line.strip().split()
            genotype_list = l[6:]
            hap_list = []
            for g1, g2 in zip(genotype_list[::2], genotype_list[1::2]):
                if g1=="N" or g2=="N" or g1=="*" or g2=="*":
                    hap_list.append("N")
                elif g1==g2:
                    hap_list.append(g1)
                else:
                    g = g1+g2
                    hap_list.append(iupac_to_nucleotides[g])
            hap_str = "".join(hap_list)
            if hap_str not in hap_dic.keys():
                hap_dic[hap_str] = [l[0]]
            else:
                hap_dic[hap_str].append(l[0])
    sorted_hap_dic = dict(sorted(hap_dic.items(), key=lambda item: len(item[1]), reverse=True))
    return sorted_hap_dic


def write_output(hap_info, output_file):
    """Write haplotype information to output file."""
    with open(output_file, 'w') as fout:
        n = 0
        for hapseq,v in hap_info.items():
            n+=1
            hap_name = "Hap"+str(n)
            fout.write(f"{hap_name}\t{len(v)}\t{hapseq}\t{','.join(v)}\n")
    print(f"Haplotype summary written to {output_file}")

def main():
    args = parse_arguments()

    # Generate haplotypes
    hap_info = generate_haplotypes(args.ped)
    # Write output
    write_output(hap_info, args.output)

if __name__ == '__main__':
    main()
