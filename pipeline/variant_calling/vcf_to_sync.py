import sys

def parse_pileup_bases(bases, ref):
    counts = {'A': 0, 'T': 0, 'C': 0, 'G': 0, 'N': 0, 'DEL': 0}
    i = 0
    while i < len(bases):
        base = bases[i]
        if base in '.,':
            counts[ref] += 1
            i += 1
        elif base == '^':
            i += 2
        elif base == '$':
            i += 1
        elif base in '+-':
            i += 1
            num = ''
            while i < len(bases) and bases[i].isdigit():
                num += bases[i]
                i += 1
            jump = int(num) if num else 0
            if base == '-': counts['DEL'] += 1
            i += jump
        elif base.upper() in 'ATCGN*':
            if base == '*': counts['DEL'] += 1
            else: counts[base.upper()] += 1
            i += 1
        else: i += 1
    return counts

def main(mpileup_file, vcf_file, sync_file):
    variant_positions = set()
    with open(vcf_file, 'r') as vcf:
        for line in vcf:
            if line.startswith('#'): continue
            f = line.strip().split('\t')
            variant_positions.add((f[0], int(f[1])))
    
    with open(mpileup_file, 'r') as mp, open(sync_file, 'w') as out:
        for line in mp:
            f = line.strip().split('\t')
            if len(f) < 6: continue
            chrom, pos, ref = f[0], int(f[1]), f[2].upper()
            if (chrom, pos) not in variant_positions: continue
            
            sync_counts = []
            num_samples = (len(f) - 3) // 3
            for i in range(num_samples):
                idx = 3 + (i * 3)
                counts = parse_pileup_bases(f[idx+1], ref)
                sync_str = f"{counts['A']}:{counts['T']}:{counts['C']}:{counts['G']}:{counts['N']}:{counts['DEL']}"
                sync_counts.append(sync_str)
            out.write(f"{chrom}\t{pos}\t{ref}\t" + "\t".join(sync_counts) + "\n")

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
