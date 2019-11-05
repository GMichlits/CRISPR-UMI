__author__ = 'georg.michlits'

input_file = open(input('enter inputfile (e.g. Ms_top50_resorted.txt): '),'r')
output_file = open(input('enter output_filename (e.g. Ms_top4.txt): '),'w')
top_X_guides = int(input('number of guides selected per gene (4): '))

genename = ''
top6 = ''
i = 0

for line in input_file:
    line_list = line.rstrip('\n').split('\t')
       #new gene write new list
    if not line_list[0] == genename:
        output_file.write(top6)
        top6 = ''
        i = 0
    genename = line_list[0]
    guideDNA = line_list[3][4:24]
    assert len(guideDNA) == 20
    if i < top_X_guides:
        top6 = top6 + line
        i = i + 1

#write last gene
output_file.write(top6)

input_file.close()
output_file.close()
print('finished')
