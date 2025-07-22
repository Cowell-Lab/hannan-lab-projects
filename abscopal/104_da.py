import airr
import csv

dirpath = '/work/data/immune/Raquib-iSABR/Raquib-abscopal/853273733430055401-242ac11d-0001-007/'
da_path = '/work/data/immune/Raquib-iSABR/Raquib-abscopal/differential_abundance/'

# patient 104
tumor_file = '24673-6.igblast.airr.tsv'
nontumor_file = '24563-5.igblast.airr.tsv'
distal_file = '26473-1.igblast.airr.tsv'

time_files = [
    "O_weeks.igblast.airr.tsv",
    "5_weeks.igblast.airr.tsv",
    "5_weeks_2.igblast.airr.tsv",
    "8_weeks.igblast.airr.tsv",
    "17_weeks.igblast.airr.tsv",
    "7_months.igblast.airr.tsv",
    "30_months.igblast.airr.tsv"
]

da1_file = '5_weeks_VS_O_weeks.differentialAbundance.tsv'
da2_file = '5_weeks_2_VS_O_weeks.differentialAbundance.tsv'
da3_file = '8_weeks_VS_O_weeks.differentialAbundance.tsv'
da4_file = '17_weeks_VS_O_weeks.differentialAbundance.tsv'
da5_file = '7_months_VS_O_weeks.differentialAbundance.tsv'
# note: the 30 month timepoint sample was named 104-D-G
da6_file = '104-D-G_VS_O_weeks.differentialAbundance.tsv'

def read_da_file(filename):
    da = []
    file = open(da_path + filename, 'r')
    row = file.readline()
    while row[0] == '#':
        row = file.readline()
    reader = csv.DictReader(file, dialect='excel-tab')
    for row in reader:
        da.append(row)
    return da
                            
def read_cdr3_dict(filename):
    cdr3 = {}
    data = airr.read_rearrangement(filename)

    for r in data:
        if not r['productive']: continue
        if cdr3.get(r['junction_aa']):
            cdr3[r['junction_aa']] = cdr3[r['junction_aa']] + r['duplicate_count']
        else:
            cdr3[r['junction_aa']] = r['duplicate_count']

    return cdr3

def venn_cdr3(na, A, nb, B):
    print(na + " total: " + str(len(A)))
    print(nb + " total: " + str(len(B)))
    X = A & B
    print("shared: " + str(len(X)))
    X = A - B
    print(na + " unique: " + str(len(X)))
    X = B - A
    print(nb + " unique: " + str(len(X)))
    print("*****")

def check_airr_file(filename, da):
    data = airr.read_rearrangement(filename)
    for row in data:
        for r in da:
            if row['sequence'] == r['sequence']:
                print(row['sequence'])
                print(row['junction_aa'], row['productive'], row['duplicate_count'], r['significance'], r['sampleA_abundance'], r['sampleB_abundance'])

def compare_airr_file(filename1, filename2, da):
    for r in da:
        found1 = False
        data = airr.read_rearrangement(filename1)
        for row in data:
            if row['sequence'] == r['sequence']:
                print(row['sequence'])
                print(row['junction_aa'], row['productive'], row['duplicate_count'], r['significance'], r['sampleA_abundance'], r['sampleB_abundance'])
                found1 = True
        found2 = False
        data = airr.read_rearrangement(filename2)
        for row in data:
            if row['sequence'] == r['sequence']:
                print(row['sequence'])
                print(row['junction_aa'], row['productive'], row['duplicate_count'], r['significance'], r['sampleA_abundance'], r['sampleB_abundance'])
                found2 = True

        if found1 and found2:
            print('FOUND IN BOTH')
        elif found1:
            print('FOUND IN 1')
        elif found2:
            print('FOUND IN 2')

da1 = read_da_file(da1_file)
da2 = read_da_file(da2_file)
da3 = read_da_file(da3_file)
da4 = read_da_file(da4_file)
da5 = read_da_file(da5_file)
da6 = read_da_file(da6_file)

print('***')
print('DA: ' + str(len(da1)))
print('TUMOR: ' + da1_file)
check_airr_file(dirpath + tumor_file, da1)
print('NON-TUMOR: ' + da1_file)
check_airr_file(dirpath + nontumor_file, da1)
print('COMPARE')
print('1: TUMOR')
print('2: NON-TUMOR')
compare_airr_file(dirpath + tumor_file, dirpath + nontumor_file, da1)

print('***')
print('DA: ' + str(len(da2)))
print('TUMOR: ' + da2_file)
check_airr_file(dirpath + tumor_file, da2)
print('NON-TUMOR: ' + da2_file)
check_airr_file(dirpath + nontumor_file, da2)
print('COMPARE')
print('1: TUMOR')
print('2: NON-TUMOR')
compare_airr_file(dirpath + tumor_file, dirpath + nontumor_file, da2)

print('***')
print('DA: ' + str(len(da3)))
print('TUMOR: ' + da3_file)
check_airr_file(dirpath + tumor_file, da3)
print('NON-TUMOR: ' + da3_file)
check_airr_file(dirpath + nontumor_file, da3)
print('COMPARE')
print('1: TUMOR')
print('2: NON-TUMOR')
compare_airr_file(dirpath + tumor_file, dirpath + nontumor_file, da3)

print('***')
print('DA: ' + str(len(da4)))
print('TUMOR: ' + da4_file)
check_airr_file(dirpath + tumor_file, da4)
print('NON-TUMOR: ' + da4_file)
check_airr_file(dirpath + nontumor_file, da4)
print('COMPARE')
print('1: TUMOR')
print('2: NON-TUMOR')
compare_airr_file(dirpath + tumor_file, dirpath + nontumor_file, da4)

print('***')
print('DA: ' + str(len(da5)))
print('TUMOR: ' + da5_file)
check_airr_file(dirpath + tumor_file, da5)
print('NON-TUMOR: ' + da5_file)
check_airr_file(dirpath + nontumor_file, da5)
print('COMPARE')
print('1: TUMOR')
print('2: NON-TUMOR')
compare_airr_file(dirpath + tumor_file, dirpath + nontumor_file, da5)

print('***')
print('DA: ' + str(len(da6)))
print('TUMOR: ' + da6_file)
check_airr_file(dirpath + tumor_file, da6)
print('NON-TUMOR: ' + da6_file)
check_airr_file(dirpath + nontumor_file, da6)
print('COMPARE')
print('1: TUMOR')
print('2: NON-TUMOR')
compare_airr_file(dirpath + tumor_file, dirpath + nontumor_file, da6)
