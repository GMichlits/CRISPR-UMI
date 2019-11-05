__author__ = 'georg.michlits'

input_file_path = (input('enter filepath input'))
input_file = open(input_file_path,'r')
out_file = open(input_file_path + str('out4mageck.txt'),'w')

i = 0
for line in input_file:
    if i == 0:
        out_file.write('sgRNA\tGene')
        column = line.rstrip('\n').split('\t')
        l_column = len(column)
        for a in range(1,len(column)):
            out_file.write('\t' + str(column[a]))
        i += 1
    else:
        column = line.rstrip('\n').split('\t')
        i += 1
        guide_string = column[0].replace('"','').rstrip(',')
        guide_list = guide_string.split(',')
        for guide in guide_list:
            gene = guide.split('_')[0]
            guide_num = guide.split('_')[1]
            guide_group = guide.split('_')[2]
            out_file.write('\n' + str(gene) + '_' + str(guide_num))
            out_file.write('\t' + str(gene))
            for a in range(1,len(column)):
                out_file.write('\t' + str(column[a]))
print('finished! ' + 'processed: ' + str(i/1000) + ' thousand guides')
