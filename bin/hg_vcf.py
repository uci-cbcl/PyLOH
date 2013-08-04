import sys

seq_string = 'ACGTacgtNn'
chrom_list = ['chrM', 'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8',
              'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15',
              'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']

def main():
    indbsnp = open(sys.argv[1])
    infasta = open(sys.argv[2])
    outfasta = open(sys.argv[3], 'w')

    snp_dict = {}

    i = 0
    for line in indbsnp:
        if line[0:3] != 'chr':
            continue

        chrom, pos, ID, ref, alt = line.strip('\n').split('\t')[0:5]

        if (len(ref) != 1) or (len(alt) != 1) or (ref not in seq_string) or (alt not in seq_string):
            continue

        seq_key = chrom + ':' + pos
        seq_value = ref + ':' + alt
        
        snp_dict[seq_key] = seq_value

        i += 1
        if i%1000000 == 0:
            print '{0} line parsed...'.format(i)

    base = 0
    for line in infasta:
        if line[0] == '>':
            chrom = line.strip('\n').strip('>')
            base = 0
            outfasta.write(line)

        else:
            seq_list = list(line.strip('\n'))
            shift = 0
            for bp in seq_list:
                idx = shift
                shift += 1
                pos = base + shift
                seq_key = chrom + ':' + str(pos)
                if snp_dict.has_key(seq_key):
                    seq_value = snp_dict[seq_key]
                    ref, alt = seq_value.split(':')
                    seq_list[idx] = alt
            line_new = "".join(seq_list) + '\n'
            outfasta.write(line_new)
            base += shift

        if base%10000000 == 0:
            print '{0}\t{1}'.format(chrom, base)

    indbsnp.close()
    infasta.close()
    outfasta.close()

if __name__ == '__main__':
    main()
