__author__ = 'georg.michlits'

##condition in shared Experiment each UMI must be present in total at least 3 times
## IMPORTANT - to avoid division by 0 the script adds pseudocount to zeros which becomes increasingly problematic if total reads of samples differ a lot.
## if total sample reads (ctrl vs exp) differ a lot try pseudocount like 0.1 or use other scripts


sample_name = input('enter input filename (e.g. MD_30kMsEtop_i0MMsg1MM.txt): ')
new_out_name = input('enter output foldername (e.g.demo ): ')
top_popped = input('enter number of clones popped (2) if no poptop (0): ')
ctrl = input('type control sample name (e.g.ctrl_200clones): ')
ctrl2 = input('type control sample #2 name (e.g.ctrl_200clones) (if no 2nd ctrl sample type:NO): ')
exp = input('type treated sample #1 name (e.g.exp_200clones): ')
exp2 = input('type treated sample #2 name (e.g.exp_200clones) (if no 2nd exp type:NO): ')
maxMM = int(input('enter max MM for removing bc-shadows (1): '))
minfold = float(input('enter min fold for removing bc-shadows (3): '))
minreads = float(input('enter minimum reads in ctrl+exp to accept clone (3): '))
pseudocount = float(input('enter pseudocount (added to 0 reads to calculate fold_change)(0.5): '))
#list_of_interesting_guide_name = input('filename of genes to look at (each line one guidename): ')

def MMmax(dna1,dna2,maxMM):
    c=0
    for i in range(0,len(dna1)):
        if not dna1[i] == dna2[i]:
            c = c + 1
            if c > maxMM:
                return False
                break
    return True

from scipy.stats import binom
import statistics
import os

input_MD = open(sample_name,'r')
current_working_directory = os.getcwd()
os.mkdir(current_working_directory + '/' + new_out_name)
#output_3_clean = open('out_3_'+sample_name + exp + 'bc' + str(maxMM)+'MM_MD_clean.txt','w')
#output_3_kicked = open('out_3_'+sample_name+ exp + 'bc' + str(maxMM)+'MM_kicked_shadows.txt','w')
intermed_4_conv = open(new_out_name + '/CONV_interm_' + new_out_name + '.txt','w')
intermed_4_CRUMI = open(new_out_name + '/CrUMI_interm_' + new_out_name + '.txt','w')
intermed_median = open(new_out_name + '/Median_interm_' + new_out_name + '.txt','w')


#create a dict/hash for
#UMIbook
# {
# guide{
#         barcode: [ctrl(int),exp(int),total(int)],
#       }
# }
#
# guidebook:{
#         ctrl: int,
#         exp: int,
#         median: float,
#         pos_clones: int,
#         total_clones: int,
#         med_of_depl: float,
#         frac_depl: float
#       }

print('reading in Masterdict')
UMIbook = {}
guidebook = {}
total_reads_exp = 0
total_reads_ctrl = 0
#every line gives number of reads for a specific guide bc combination and specific exp

for line in input_MD:
    element = line.rstrip('\n').split('\t')
    guide = element[0].rstrip(',')
    exp_name = element[1]
    barcode = element[2]
    count = int(element[3])

# here i limit the script to the experiments mentioned in the input to be analysed
    #### for classical analysis guidebook[guide] = Information
    #### for CRISPR DigDeep analysis UMIbook[guide][barcode] = (ctr count, exp count, sum of both count)
#  only later i go through the data to remove shadows etc.

    if exp_name == ctrl or exp_name == ctrl2:
        total_reads_ctrl += count
        if not guide in guidebook:
            guidebook[guide] = {'ctrl':0, 'exp':0, 'median':0, 'pos_clones':0,
                                'med_of_depl':0, 'frac_depl':0}
        guidebook[guide]['ctrl'] += count
        if not guide in UMIbook:
            UMIbook[guide] = {}
        if not barcode in UMIbook[guide]:
            UMIbook[guide][barcode] = [0,0,0]   #[0] ctrl counts, [1], exp counts, [2] total counts
        #0 is the index for ctrl counts
        UMIbook[guide][barcode][0] = count
        #3 is the index for total reads
        UMIbook[guide][barcode][2] += count
    if exp_name == exp or exp_name == exp2:
        total_reads_exp += count
        if not guide in guidebook:
            guidebook[guide] = {'ctrl':0, 'exp':0, 'median':0, 'pos_clones':0,
                                'med_of_depl':0, 'frac_depl':0}
        guidebook[guide]['exp'] += count
        if not guide in UMIbook:
            UMIbook[guide] = {}
        if not barcode in UMIbook[guide]:
            UMIbook[guide][barcode] = [0,0,0]
        #1 is the index for exp counts
        UMIbook[guide][barcode][1] = count
        #3 is the index for total reads
        UMIbook[guide][barcode][2] += count

print('reading in MD finished! - start shadow-kicking and calculations')
print('total_reads_ctrl =' + str(total_reads_ctrl))
print('total_reads_exp =' + str(total_reads_exp))

#insert poptop algorithm very simple it takes a list of guides and removes the top X (2) clones

poptop_report = open(new_out_name + '/Report_poptop' + new_out_name + '.txt','w')
#script pops the top X clones for all guides
poptop_report.write('sgRNA\tbarcode_pop\tctrl_reads\texp_reads\n')
x=0
for guide in sorted(UMIbook):
    a_pop = 0
    for barcode_item in sorted(UMIbook[guide].items(), key=lambda t: (t[1][0] + t[1][1]), reverse=True):
        if a_pop <= (int(top_popped)-1):
            popped_ctr_reads = barcode_item[1][0]
            total_reads_ctrl = total_reads_ctrl - popped_ctr_reads
            popped_exp_reads = barcode_item[1][1]
            total_reads_exp = total_reads_exp - popped_exp_reads
            UMIbook[guide].pop(barcode_item[0])
            poptop_report.write(guide + '\t' + str(barcode_item[0]) + '\t' + str(barcode_item[1][0]) + '\t' + str(barcode_item[1][1]) + '\n')
            x += 1
        a_pop += 1
print('popped ' + str(x) + ' clones.')
print('total_reads_ctrl =' + str(total_reads_ctrl))
print('total_reads_exp =' + str(total_reads_exp))


# generate a new UMIbook that is free of barcode_shadows (i use barcode as a synonym for UMI - the 10nt sequences identifying clones)
#UMIbook_clean
# {
# guide{
#         barcode: [ctrl(int),exp(int),total_ctrl(int),totalexp(),pctrl/totalctr,rpm_ctrl,rpm_exp,pdepl,penri,fold_c,p_vulcano],
#       }
# }

# generate a new (dict/hash) Medianhash to find fold_c cutoff
trash_count = 0
trash_count_reads = 0
total_bc = 0
total_bc_reads = 0
UMIbook_clean = {}
Medianhash = {}
list_num_clones = []
list_avreadsclones = []
for guide in sorted(UMIbook):   #for every guide....
    ref_barcode_list = []
    clone_fold_list = []
    guide_reads_list = []
    for item in sorted(UMIbook[guide].items(), key=lambda t: (t[1][1]+t[1][0]), reverse=True):    #sorted to begin with clone with most reads
        barcode = item[0]
        total_reads_barcode = float(item[1][1]) + float(item[1][0])
        total_bc += 1
        total_bc_reads += total_reads_barcode
        shadow_found = 'no'  #set to 'no before searching
        if total_reads_barcode < minreads:                  #input variable minreads, minimium read a clone needs to have to accept it in analysis
            shadow_found = 'yes'
            trash_count += 1
            trash_count_reads += total_reads_barcode
        else:
            ref_barcode_list.append((barcode,total_reads_barcode))
            for ref_barcode in ref_barcode_list:
                ref_bc_seq = ref_barcode[0]
                ref_bc_count = ref_barcode[1]
                if (MMmax(barcode,ref_bc_seq,maxMM) and total_reads_barcode * minfold <= ref_bc_count):
                    shadow_found = 'yes'
                    #(guide + 'shadow eliminated')
                    trash_count += 1
                    trash_count_reads += total_reads_barcode
                    #here one could read in a shadow book to then later print out bc_shadows sorted out
                    break
        if shadow_found == 'no':
            temp_list = UMIbook[guide][barcode] #temp_list is now [0] ctrl counts, [1], exp counts, [2] total counts (exp+ctrl)
            ctrl_reads = temp_list[0]
            exp_reads = temp_list[1]
            temp_list[2] = total_reads_ctrl #2 now: [0] ctrl counts, [1], exp counts, [2] total counts ctrl whole exp,
            temp_list.append(total_reads_exp) #3 now: [0] ctrl counts, [1], exp counts, [2] total counts ctrl whole exp, [3] total count exp [whole exp]
            p_ctr_totctr = ctrl_reads / total_reads_ctrl
            if ctrl_reads == 0:
                p_ctr_totctr = pseudocount / total_reads_ctrl
            temp_list.append(p_ctr_totctr) #4
            rpm_ctrl = ctrl_reads / total_reads_ctrl
            if ctrl_reads == 0:
                rpm_ctrl = pseudocount / total_reads_ctrl
            temp_list.append(rpm_ctrl) #5
            rpm_exp = exp_reads / total_reads_exp
            if exp_reads == 0:
                rpm_exp = pseudocount / total_reads_exp
            temp_list.append(rpm_exp) #6
            p_value_depl = (binom.cdf(exp_reads,total_reads_exp,p_ctr_totctr))
            p_value_enri = (binom.sf(exp_reads,total_reads_exp,p_ctr_totctr))
            temp_list.append(p_value_depl) #7
            temp_list.append(p_value_enri) #8
            fold_c = rpm_exp / rpm_ctrl
            temp_list.append(fold_c) #9
            if fold_c <= 1:
                p_vulcano = p_value_depl
            else:
                p_vulcano = p_value_enri
            temp_list.append(p_vulcano) #10
            if not guide in UMIbook_clean:
                UMIbook_clean[guide] = {}
            UMIbook_clean[guide][barcode] = temp_list

            #write fold-c in a list for that guide #within that guide
            clone_fold_list.append(fold_c)
            guide_reads_list.append(total_reads_barcode/2) #because it represents 2 conditions exp and ctrl and i want to know average reads per condtion
    number_of_clones = len(clone_fold_list)
    if not number_of_clones == 0:
        median_fold_c = statistics.median(clone_fold_list)
        average_reads_perclone = statistics.mean(guide_reads_list)
        list_num_clones.append(number_of_clones)
        list_avreadsclones.append(average_reads_perclone)
    else:
        list_num_clones.append(0)
        list_avreadsclones.append(0)

    if number_of_clones >= 5:
        Medianhash[guide] = [median_fold_c,number_of_clones,average_reads_perclone]
av_clones = statistics.mean(list_num_clones)
av_reads_perUMI = statistics.mean(list_avreadsclones)
print('finished reading in MD')
print('total UMI analysed: ' + '\t' + str(total_bc))
print('total UMI trashed (shadows): ' + '\t' + str(trash_count/total_bc*100) + '%\t' + 'left with' + str(total_bc-trash_count) + 'clones')
print('total reads: ' + '\t' + str(total_bc_reads) + '\t' + 'left with ' + str((total_bc_reads-trash_count_reads)/total_bc_reads*100) + '% of reads')
print('average clones per guide: ' + '\t' + str(av_clones) + '\t' + 'average reads per clone: ' + '\t' + str(av_reads_perUMI))
clones_in_analysis = total_bc-trash_count

print('writing result files')
# write guide score by median
i = 0

intermed_median.write('sgRNA\tmedian_(treated/ctrl)\tnumber_of_clones\taverge_reads_per_clone\n')

for item in sorted(Medianhash.items(), key=lambda t: t[1][0]): #sorted by median_fold_c
    guide = item[0]
    median_fold_c = item[1][0]
    number_of_clones = item[1][1]
    average_reads_perclone = item[1][2]
    intermed_median.write(guide + '\t' + str(median_fold_c) + '\t' + str(number_of_clones) + '\t' + str(average_reads_perclone) + '\n')
intermed_median.close()

#now i run a small script to fix some issues with guide names (the same guides in different subpools have different names, at the same time i do not consider genes from the predicted off-target library subpool)

intermed_median = open(new_out_name + '/Median_interm_' + new_out_name + '.txt','r')
out_median_unique = open(new_out_name + '/CrUMI_median_' + new_out_name + '.txt','w')

n = 0
for line in intermed_median:
    if n == 0:
        out_median_unique.write(line)
    column = line.rstrip('\n').split('\t')
    name = column[0]
    del column[0]
    name_list = name.split(',')
    new_names = []
    for element in name_list:
        if (not 'OfT' in element) and '_' in element:
            guide = element.split('_')[0]+'_'+element.split('_')[1]
            if guide not in new_names:
                new_names.append(guide)
    for shortname in new_names:
        out_median_unique.write(shortname)
        for thing in column:
            out_median_unique.write('\t' + thing)
        out_median_unique.write('\n')
    n += 1
print('Median_file finished.' + '\t' + 'processed: ' + str(n) + ' guides')
out_median_unique.close()


#write details on clones for CRISPR-UMI Analysis

intermed_4_CRUMI.write('sgRNA\tbarcode\tctrl_reads\texp_reads\ttotal_ctrl\ttotal_exp\tp_binom\trpm_ctrl\trpm_exp\tp_depl\tp_enri\tfold_c\tp_vulcano\n')
#ctrl(int),exp(int),total_ctrl(int),totalexp(),pctrl/totalctr,rpm_ctrl,rpm_exp,pdepl,penri,fold_c,p_vulcano
for guide in sorted(UMIbook_clean):
    for bc in UMIbook_clean[guide]:
        intermed_4_CRUMI.write(guide + '\t' + bc)
        for element in UMIbook_clean[guide][bc]:
            intermed_4_CRUMI.write('\t' + str(element))
        intermed_4_CRUMI.write('\n')
intermed_4_CRUMI.close()

#now i run a small script to fix some issues with guide names (the same guides in different subpools have different names, at the same time i do not consider genes from the predicted off-target library subpool)

intermed_4_CRUMI = open(new_out_name + '/CrUMI_interm_' + new_out_name + '.txt','r')
out_Crumi_unique = open(new_out_name + '/CrUMI_unique_' + new_out_name + '.txt','w')
n = 0
for line in intermed_4_CRUMI:
    if n == 0:
        out_Crumi_unique.write(line)
    column = line.rstrip('\n').split('\t')
    name = column[0]
    del column[0]
    name_list = name.split(',')
    new_names = []
    for element in name_list:
        if (not 'OfT' in element) and '_' in element:
            guide = element.split('_')[0]+'_'+element.split('_')[1]
            if guide not in new_names:
                new_names.append(guide)
    for shortname in new_names:
        out_Crumi_unique.write(shortname)
        for thing in column:
            out_Crumi_unique.write('\t' + thing)
        out_Crumi_unique.write('\n')
    n += 1
out_Crumi_unique.close()


#now a short script that generates an input file for MAGeCK

out_Crumi_unique = open(new_out_name + '/CrUMI_unique_' + new_out_name + '.txt','r')
out_Crumi_mageck = open(new_out_name + '/CrUMI_4mageck_' + new_out_name + '.txt','w')
i = 0
for line in out_Crumi_unique:
    column = line.rstrip('\n').split('\t')
    if i == 0:
        out_Crumi_mageck.write('sgRNA\tGene\tctr\texp\n')
    else:
        column = line.rstrip('\n').split('\t')
        guide = column[0]
        clone = guide + '_' + column[1]
        ctrl = column[2]
        exp = column[3]
        if not (int(ctrl) + int(exp)) == 0:
            out_Crumi_mageck.write(clone + '\t' + guide + '\t' + ctrl + '\t' + exp + '\n')#'sgRNA\tGene\tc8_cl\te8_cl\n')
    i+=1

print('CrUMI_file finished ' + 'processed: ' + str(i) + ' clones')


# write out conventional
#UMIbook
# {
# guide{
#         barcode: [ctrl(int),exp(int),total(int)],

intermed_4_conv.write('sgRNA\tctr\texp\tpdepl\tpenri\tfold_c\tp_vulcano\n')
for guide in UMIbook:
    intermed_4_conv.write(guide)
    contr_exp_list = [0,0]
    for bc in UMIbook[guide]:
        contr_exp_list[0] += UMIbook[guide][bc][0]
        contr_exp_list[1] += UMIbook[guide][bc][1]
    ctrl = contr_exp_list[0]
    rpm_ctrl = ctrl / total_reads_ctrl*1000000
    if ctrl == 0:
        rpm_ctrl = pseudocount / total_reads_ctrl*1000000
    exp = contr_exp_list[1]
    rpm_exp = exp / total_reads_exp*1000000
    if exp == 0:
        rpm_exp = pseudocount / total_reads_exp*1000000
    p_ctr_totctr = ctrl / total_reads_ctrl
    if ctrl == 0:
        p_ctr_totctr = pseudocount / total_reads_ctrl
    p_value_depl = (binom.cdf(exp,total_reads_exp,p_ctr_totctr))
    p_value_enri = (binom.sf(exp-1,total_reads_exp,p_ctr_totctr))
    fold_c = rpm_exp / rpm_ctrl
    if fold_c <= 1:
        p_vulcano = p_value_depl
    else:
        p_vulcano = p_value_enri
    intermed_4_conv.write(guide + '\t' + str(ctrl) + '\t' + str(exp) + '\t' + str(p_value_depl) + '\t' + str(p_value_enri) +'\t' + str(fold_c) +'\t' + str(p_vulcano)+ '\n')
intermed_4_conv.close()

#now i run a small script to fix some issues with guide names (the same guides in different subpools have different names, at the same time i do not consider genes from the predicted off-target library subpool)

intermed_4_conv = open(new_out_name + '/CONV_interm_' + new_out_name + '.txt','r')
out_conv_unique = open(new_out_name + '/CONV_unique_' + new_out_name + '.txt','w')
n = 0
for line in intermed_4_conv:
    if n == 0:
        out_conv_unique.write(line)
    column = line.rstrip('\n').split('\t')
    name = column[0]
    del column[0]
    name_list = name.split(',')
    new_names = []
    for element in name_list:
        if (not 'OfT' in element) and '_' in element:
            guide = element.split('_')[0]+'_'+element.split('_')[1]
            if guide not in new_names:
                new_names.append(guide)
    for shortname in new_names:
        out_conv_unique.write(shortname)
        for thing in column:
            out_conv_unique.write('\t' + thing)
        out_conv_unique.write('\n')
    n += 1

out_conv_unique.close()

#now a short script that generates an input file for MAGeCK

out_conv_unique = open(new_out_name + '/CONV_unique_' + new_out_name + '.txt','r')
out_conv_Mageck = open(new_out_name + '/CONV_4mageck_' + new_out_name + '.txt','w')

i = 0
for line in out_conv_unique:
    column = line.rstrip('\n').split('\t')
    if i == 0:
        out_conv_Mageck.write('sgRNA\tGene\tctr\texp\n')
    else:
        sgRNA = column[0].replace('"','')
        Gene = column[0].replace('"','').split('_')[0]
        c200 = column[1]
        e200 = column[2]
        out_conv_Mageck.write(sgRNA + '\t' + Gene + '\t' + c200 + '\t' + e200 + '\n')
    i += 1
print('CONV_file finished!' + '\t' 'processed: ' + str(i) + ' guides')
print('finished')
