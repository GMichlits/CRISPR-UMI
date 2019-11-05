__author__ = 'georg.michlits'

import re

input_filename = input('enter inputfilename (e.g. mergeData.2MM.new.txt): ')
output_filename = input('enter outputfilename (e.g. Ms_top50_resorted.txt): ')
list10genomewide = open(input_filename, 'r')
top10termback = open(output_filename,'w')

genename = ''
list_back = ''
Terminatorless = ''
Terminator_in = ''
so_l = 0
sortout_lis = ''
t = 0

for line in list10genomewide:
    line_list = line.rstrip('\n').split('\t')
    #new gene write new list
    if not line_list[0] == genename:
        # print previous gene list terminator_in last
        list_back = Terminatorless + Terminator_in
        top10termback.write(list_back)
        list_back = ''
        Terminatorless = ''
        Terminator_in = ''
    genename = line_list[0]
    guideDNA = line_list[3][4:24]
    guidestart = line_list[3][4:9]
    assert len(guideDNA) == 20
    if guidestart == 'TCTTC':
        sortout_lis = sortout_lis + line
        so_l = so_l + 1
    elif not re.search(r'T{4,20}',guideDNA):
        Terminatorless = Terminatorless + line
    else:
        Terminator_in = Terminator_in + line
        t = t + 1

#list terminator_in last
list_back = Terminatorless + Terminator_in
top10termback.write(list_back)

list10genomewide.close()
top10termback.close()
print('finished')