__author__ = 'georg.michlits'

filename = input('type filename (.csv)')
file = open(filename,'r')
file_out = open(filename + '_out','w')

import statistics

genedict = {}
scorelist = []
i=0
for line in file:
    if not line.startswith('sgRNA'):
        i += 1
        column = line.rstrip('\n').split(',')
        guide_name = column[0]
        genename = guide_name.split('_')[0]
        if not genename in genedict:
            genedict[genename] = 0
        genedict[genename] += 1
        for element in genedict:
            scorelist.append(genedict[element])
        mean_scorelist = statistics.mean(scorelist)
        file_out.write(guide_name + ',' + column[1] + ',' + str(i) + ',' + str(len(scorelist)) + ',' + str(mean_scorelist) + '\n')
        scorelist = []
file.close()
file_out.close()


