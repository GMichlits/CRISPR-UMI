__author__ = 'georg.michlits'

def MMmax(dna1,dna2,maxMM):
    c=0
    for i in range (0,len(dna1)):
        if not dna1[i] == dna2[i]:
            c = c + 1
            if c > maxMM:
                return False
                break
    return True

sample_file_name = input('enter filename Masterdict (e.g. iPSC_MD): ')
maxMM = int(input('enter max MM for removing bc-shadows (1): '))
minfold = float(input('enter min fold for removing bc-shadows (10): '))
minreads = float(input('enter minimum reads in ctrl+exp to accept clone (10): '))
new_name = input('enter new sample name (e.g. iPSC_1MM10X10min): ')
input_3_out2_Masterdict = open(sample_file_name,'r')
output_3_clean = open('out_3_'+ new_name + '_MDclean.txt','w')
output_4_class = open('out_4_'+ new_name + '_conv.txt','w')
output_3_kicked = open('out_3_'+ new_name + '_MDkicked.txt','w')
indexfilename = input('enter indexfilename (e.g. ind_file.csv): ')
Indexfile = open(indexfilename,'r')

#readout_class = {'gDNAname:{index:(Incidence,Abundance)}}
#readout_bc = {'gDNAname:{index:{barcode:count}}}
readout_shadows = {}
readout_class = {}
readout_bc = {}
count_sum = 0
barcode_list = []
bc_list_classreadout = []
sh_pro = 0
sh_el = 0
incidence = 0

print('reading in Masterdict')
#every line gives number of reads for a specific guide bc combination and specific exp
for line in input_3_out2_Masterdict:
    sh_pro = sh_pro + 1
    element = line.rstrip('\n').split('\t')
    if sh_pro == 1:
        gDNAname = element[0]
        exp_name = element[1]
        barcode = element[2]
        count = int(element[3])
    gDNAname_n = element[0]
    exp_name_n = element[1]
    barcode_n = element[2]
    count_n = int(element[3])
    if gDNAname_n+exp_name_n == gDNAname+exp_name: #same guide same exp
        last_line_new_gene = 'no'
        shadow_found = 'no'     #resets to 'no' - before checking
        count_sum = count_sum + count_n
        for bc_c_tuple in barcode_list:        #now check reference list
            ref_barcode = bc_c_tuple[0]
            ref_count = bc_c_tuple[1]

            if count_n < minreads or (MMmax(barcode_n,ref_barcode,maxMM) and count_n * minfold <= ref_count):
                shadow_found = 'yes'
                barcode_list.append((barcode_n,count_n))
                output_3_kicked.write(line)
                sh_el = sh_el + 1
                if not gDNAname_n in readout_shadows:
                    readout_shadows[gDNAname_n] = {}
                if not barcode_n in readout_shadows[gDNAname_n]:
                    readout_shadows[gDNAname_n][barcode_n] = {}
                readout_shadows[gDNAname_n][barcode_n][exp_name_n] = count_n
                #print('shadows eliminated: ' + str(sh_el) + '\t' + 'lines processed:' + str(sh_pro) + '\t' + str(count_n))
                break
        if shadow_found == 'no':
            incidence = incidence + 1
            barcode_list.append((barcode_n,count_n))
            output_3_clean.write(line)
            if not gDNAname_n in readout_bc:
                readout_bc[gDNAname_n] = {}
            if not barcode_n in readout_bc[gDNAname_n]:
                readout_bc[gDNAname_n][barcode_n] = {}
            readout_bc[gDNAname_n][barcode_n][exp_name_n] = count_n

    else:       # new exp index - or new guide
        #write data for old genename_exp
        # add total counts print('len' + str(len(barcode_list)))
        output_4_class.write(gDNAname + '\t' + exp_name + '\t' + str(incidence) + '\t' + str(count_sum) + '\n')
        if not gDNAname in readout_class:
            readout_class[gDNAname] = {}
        readout_class[gDNAname][exp_name] = (incidence,count_sum)
        #write data for new line "always a true bc"
        if not gDNAname_n in readout_bc:
            readout_bc[gDNAname_n] = {}
        if not barcode_n in readout_bc[gDNAname_n]:
            readout_bc[gDNAname_n][barcode_n] = {}
        readout_bc[gDNAname_n][barcode_n][exp_name_n] = count_n
        #reset conditions for new genename_exp
        if count_n >= minreads:
            incidence = 1
        else:
            incidence = 0
        count_sum = count_n
        barcode_list = []
        barcode_list.append((barcode_n,count_n))
        output_3_clean.write(line)
        gDNAname = gDNAname_n
        exp_name = exp_name_n
        barcode = barcode_n
        count = count_n
#WRITE OUT LAST LINE
output_4_class.write(gDNAname + '\t' + exp_name + '\t' + str(incidence) + '\t' + str(count_sum) + '\n')
if not gDNAname in readout_class:
    readout_class[gDNAname] = {}
readout_class[gDNAname][exp_name] = (incidence,count_sum)
output_4_class.close()
output_3_kicked.close()
output_3_clean.close()
print('reading in done')
print('finding bc-shadows done')
print('shadows eliminated: ' + str(sh_el) + '\t' + 'lines processed:' + str(sh_pro))

explist = []
for line in Indexfile:
    element = line.rstrip('\n').split(',')
    Exp_name = element[1]
    if not Exp_name in explist:
        explist.append(Exp_name)
print('Index dictionaries done')
Indexfile.close()

#write output 4 (MD style into table style)
output_4_class = open('out_4_'+ new_name + '_conv.txt','r')
output_5_class = open('out_5_'+ new_name + '_IncAbu.txt','w')

#read in dict for 4class
classdict ={}
i = 0
for line in output_4_class:
    i += 1
    element = line.rstrip('\n').split('\t')
    gDNAname = element[0]
    exp_name = element[1]
    Incidence = element[2]
    Abundance = element[3]
    if not gDNAname in classdict:
        classdict[gDNAname] = {}
    classdict[gDNAname][exp_name] = (Incidence,Abundance)
#write output_5_class
print('wrote ' + str(i) + ' guides in classical output file')

#header
output_5_class.write('guidename')
v = 'In'
for exp in explist:
    output_5_class.write('\t' + exp + '_' + v)
v = 'Ab'
for exp in explist:
    output_5_class.write('\t' + exp + '_' + v)
output_5_class.write('\n')

#table
for gDNAname in sorted(classdict):
    output_5_class.write(gDNAname)
    #give index for Incidence
    index_var = 0
    for exp in explist:
        if exp in classdict[gDNAname]:
            output_5_class.write('\t' + str(classdict[gDNAname][exp][index_var]))
        else:
            output_5_class.write('\t' + str(0))
    #give index for Abundance
    index_var = 1
    for exp in explist:
        if exp in classdict[gDNAname]:
            output_5_class.write('\t' + str(classdict[gDNAname][exp][index_var]))
        else:
            output_5_class.write('\t' + str(0))
    output_5_class.write('\n')


print('writing read tables')
output_5_table = open('out_6_' + new_name + '_CrUMI_' + '.txt','w')
#write header
output_5_table.write('guide\tbarcode')
for exp_name in explist:
    output_5_table.write('\t' + exp_name)
output_5_table.write('\n')

print('write output_5 per clone analysis')

for gDNAname in sorted(readout_bc):
    for barcode in readout_bc[gDNAname]:
        output_5_table.write(gDNAname + '\t' + barcode)
        i = 0
        for exp_name in explist:
            i = i + 1
            #print (i)            that is just spitting out 0 results if guides barcodes experiments are not presen
            if not exp_name in readout_bc[gDNAname][barcode]:
                if not gDNAname in readout_shadows:
                    readout_shadows[gDNAname] = {}
                if not barcode in readout_shadows[gDNAname]:
                    readout_shadows[gDNAname][barcode] = {exp_name: 0}
                if not exp_name in readout_shadows[gDNAname][barcode]:
                    c = 0
                else:
                    c = int(readout_shadows[gDNAname][barcode][exp_name])
            else:
                c = int(readout_bc[gDNAname][barcode][exp_name])
            output_5_table.write('\t' + str(c))
        output_5_table.write('\n')
print('finished')
