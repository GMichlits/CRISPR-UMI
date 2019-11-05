__author__ = 'georg.michlits'

input_file = open(input('enter input filename (e.g. Ms_druggable_with_negctrl.txt): '),'r')
output_file = open(input('enter output filename (e.g. Ms_druggable_order.txt): '),'w')
upstream = input('enter upstream sequence (e.g. GATTACATGGTCAGACGAAGACgaCACCG): ')
downstream = input('enter downstream Sequence (e.g. GTTTctGTCTTCTTACCACACCAGTCGA): ')
subpool_name = input('enter subpool name (e.g. Ms_Druggable): ')

guide_num = 0
total_guides = 0
genename = 'anything to start not a genename'
genecount = 0

for line in input_file:
    element = line.rstrip('\n').split('\t')
    genename_n = element[0]
    if genename_n == genename:
        guide_num = guide_num + 1
        total_guides = total_guides + 1
    else:
        guide_num = 1
        total_guides = total_guides + 1
        genecount = genecount + 1
    genename = genename_n
    pos = element[1]
    groupsize = 'n'
    guide_20nt = element[3][4:24]
    uniquename = subpool_name + '_' + genename_n + '_' + str(guide_num)
    output_file.write(genename_n + '_' + str(guide_num) + '\t' + pos + '\t' + subpool_name + '\t' + groupsize + '\t' + upstream + '\t' + guide_20nt + '\t' + downstream + '\t' + uniquename + '\t' + upstream.upper() + guide_20nt.upper() + downstream.upper() + '\n')

input_file.close()
output_file.close()

print('total genes: ' + str(genecount) + ' including ctrl guides')
print('total guides: ' + str(total_guides) + ' including ctrl guides')

