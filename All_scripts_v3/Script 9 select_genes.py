__author__ = 'georg.michlits'

target_genes = open(input('enter filename of target gene_list (e.g. druggable_genes_toronto.ms.txt): '), 'r')
allgenes_top4 = open(input('input filename (e.g. Ms_top4.txt): '), 'r')
output_target_guides = open(input('enter output filename (e.g. Ms_druggable.txt): '), 'w')
x = int(input('number of guides that go into screen. (4): '))

i = 1
all_guides = {}
genename = 'startname'
for line in allgenes_top4:
    line_list = line.rstrip('\n').split('\t')
    genename_n = line_list[0]
    if genename_n == genename:
        i = i + 1
        guide = line_list[0] + '_' + str(i)
        all_guides[guide] = line
    else:
        i = 1
        guide = line_list[0] + '_' + str(i)
        genename = genename_n
        all_guides[guide] = line
allgenes_top4.close()

target_guide_list = []
#read in gene_list
for line in target_genes:
    for i in range(1,x+1):  #number of guides per gene come in here again
        element = line.rstrip('\n').split('\t')
        guide = element[0] + '_' + str(i)
        target_guide_list.append(guide)
target_genes.close()

not_found = 0
found = 0
for guide in target_guide_list:
    if guide in all_guides:
        output_target_guides.write(all_guides[guide])
        found = found + 1
    else:
        not_found = not_found + 1

print('guides available: ' + str(found))
print('guides not available: ' + str(not_found))
print('finished')
