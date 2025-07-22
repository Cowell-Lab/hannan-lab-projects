import airr
import csv

dirpath = '/work/data/immune/Raquib-iSABR/Raquib-abscopal/853273733430055401-242ac11d-0001-007/'
da_path = '/work/data/immune/Raquib-iSABR/Raquib-abscopal/differential_abundance/'

# patient 107
tumor_file = 'Thr_107_31626-3.igblast.airr.tsv'
nontumor_file = 'Thr_107_31626_1.igblast.airr.tsv'
distal_file = 'Thr_107_31627-1.igblast.airr.tsv'

time_files = [
    '107_0_week.igblast.airr.tsv',
    '107_3_month.igblast.airr.tsv'
]

da_file = '107_3_month_VS_107_0_week.differentialAbundance.tsv'

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
    for r in da:
        data = airr.read_rearrangement(filename)
        for row in data:
            if row['sequence'] == r['sequence']:
                print(row['sequence'])
                print(row['junction_aa'], row['productive'], row['duplicate_count'], r['significance'])

da_107 = read_da_file(da_file)

print('TUMOR')
check_airr_file(dirpath + tumor_file, da_107)

print('NON-TUMOR')
check_airr_file(dirpath + nontumor_file, da_107)
