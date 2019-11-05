__author__ = 'georg.michlits'

filename_1 = input('enter filename_2 (path): ')
filename_2 = input('enter filename_1 (path): ')
file_merge_name = input('enter merged file_filename: ')
file1 = open(filename_1,'r')
file2 = open(filename_2,'r')
file_merge = open(file_merge_name,'w')

print('reading in file1')
for line in file1:
    file_merge.write(line)
print('finished')
print('reading in file2')
for line in file2:
    file_merge.write(line)
file1.close()
file2.close()
file_merge.close()
print('finished')
