__author__ = 'georg.michlits'

clone_list_file = open(input('enter filename e.g. CrUMI_interm.txt'),'r')
list_of_interesting_guide_name = input('filename of genes to look at (each line one guidename): ')
look_closer_file = open(list_of_interesting_guide_name,'r')
output = open('clones_'+list_of_interesting_guide_name + '.txt','w')

selection = {}
for line in look_closer_file:
    guide = line.rstrip('\n')
    print(guide)
    selection[guide] = 0
i=0
for line in clone_list_file:
    if i == 0:
        e=0
    else:
        column = line.rstrip('\n').split('\t')
        if '_' in column[0] or '':
            guidename = column[0].split('_')[0]+'_'+column[0].split('_')[1]
            #print(guidename)
        if column[0] == '':
            guidename = guidename_old
        guidename_old = guidename
        if guidename in selection:
            output.write(guidename + '\t' + line)
                #print(guidename)
    i += 1
print('finished')




