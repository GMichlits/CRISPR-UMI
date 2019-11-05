__author__ = 'georg.michlits'

# easy script that just adds together the column Etoposide 8 days for all guides of the same gene and gives p values next to each other so that
# they can be worked with easily in excel.


filename = input('type filename: ')
infile = open(filename,'r')
outfile = open(filename + 'out.tsv','w')

genename = 'test'
depl_list = []
pvalue_list = []

i = -1
for line in infile:
    element = line.rstrip('\n').split('\t')
    genename_n = element[0].split('_')[0]
    depl = float(element[1])
    pvalue = float(element[2])
    i = i + 1
    if genename_n == genename:
        depl_list.append(depl)
        pvalue_list.append(pvalue)
    else:
        #print(genename + '\t' +str(i))
        outfile.write(genename + '\t')
        len_i = i
        if not i >= 5:
            i = i + 1
            depl_list.append(1)
            pvalue_list.append(1)
        if not i >= 5:
            i = i +1
            depl_list.append(1)
            pvalue_list.append(1)
        if not i >= 5:
            i = i +1
            depl_list.append(1)
            pvalue_list.append(1)
        if not i >= 5:
            i = i +1
            depl_list.append(1)
            pvalue_list.append(1)
        if not i >= 5:
            i = i +1
            depl_list.append(1)
            pvalue_list.append(1)
        for e in sorted(pvalue_list):
            outfile.write(str(depl_list[pvalue_list.index(e)]) + '\t')
        for e in sorted(pvalue_list):
            outfile.write(str(e) + '\t')
        outfile.write(str(len_i) + '\n')
        i = 0
        depl_list = []
        pvalue_list = []
        depl_list.append(depl)
        pvalue_list.append(pvalue)
        genename = genename_n
outfile.write(genename + '\t')
i = i + 1
#print(genename + '\t' +str(i))
if not i >= 5:
    i += 1
    depl_list.append(1)
    pvalue_list.append(1)
if not i >= 5:
    i += 1
    depl_list.append(1)
    pvalue_list.append(1)
if not i >= 5:
    i += 1
    depl_list.append(1)
    pvalue_list.append(1)
if not i >= 5:
    i += 1
    depl_list.append(1)
    pvalue_list.append(1)
if not i >= 5:
    i += 1
    depl_list.append(1)
    pvalue_list.append(1)
for e in sorted(pvalue_list):
    outfile.write(str(depl_list[pvalue_list.index(e)]) + '\t')
for e in sorted(pvalue_list):
    outfile.write(str(e) + '\t')
outfile.write(str(len_i) + '\n')
