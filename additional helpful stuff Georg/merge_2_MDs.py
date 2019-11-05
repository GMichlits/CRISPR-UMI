__author__ = 'georg.michlits'

file1 = open(input('MD1: ') ,'r')
file2 = open(input('MD2: ') ,'r')
fileout = open(input('MD_out: ') ,'w')

readout = {}
for line in file1:
    element = line.rstrip('\n').split('\t')
    if len(element) == 4:
        gDNAname = element[0]
        exp_name = element[1]
        bc = element[2]
        c = int(element[3])
        if not gDNAname in readout:
            readout[gDNAname] = {}
        if not exp_name in readout[gDNAname]:
            readout[gDNAname][exp_name] = {}
        readout[gDNAname][exp_name][bc] = c

print('finished MD 1')

for line in file2:
    element = line.rstrip('\n').split('\t')
    if len(element) == 4:
        gDNAname = element[0]
        exp_name = element[1]
        bc = element[2]
        c = int(element[3])
        if not gDNAname in readout:
            readout[gDNAname] = {}
        if not exp_name in readout[gDNAname]:
            readout[gDNAname][exp_name] = {}
        if bc in readout[gDNAname][exp_name]:
            readout[gDNAname][exp_name][bc] += c
        else:
            readout[gDNAname][exp_name][bc] = c
print('finished MD 2')

print('writing Masterdict')
for gDNAname in sorted(readout):
    for exp_name in sorted(readout[gDNAname]):
        for i in sorted(readout[gDNAname][exp_name].items(), key=lambda t: t[1], reverse=True):
            barcode = i[0]
            c = i[1]
            fileout.write(gDNAname + '\t' + exp_name + '\t' + barcode + '\t' + str(c) + '\n')

print('writing Masterdict finished')
fileout.close()