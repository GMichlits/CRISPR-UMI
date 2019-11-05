__author__ = 'georg.michlits'

file1 = open(input('enter filename 1 to compare guidename in column[0] (e.g. CrUMI.gene.summary.txt_): '),'r')
file2 = open(input('enter filename 2 guidename in column[0] (e.g. CrUMI_median_pop2_200.txt) '),'r')
combo_file = open(input('enter filename of combined file (e.g. CrUMI_guides): '),'w')

n = 0
DICT = {}
for line in file1:
    if n == 0:
        header = line.rstrip('\n')
        #dosomething
    else:
        column = line.rstrip('\n').split('\t')
        name = column[0]
        DICT[name] = line.rstrip('\n')
    n += 1

k = 0
m = 0
for line in file2:
    if k == 0:
        header_2 = line
        combo_file.write(header + '\t' + header_2)
    else:
        column = line.rstrip('\n').split('\t')
        name = column[0]
        if name in DICT:
            combo_file.write(DICT[name] + '\t' + line)
            m+=1
    k += 1
print('lines in file1: \t' + str(n))
print('lines in file2: \t' + str(k))
print('matching guides: \t' + str(m))
print('finished')