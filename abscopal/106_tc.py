import airr

dirpath = '/work/data/immune/Raquib-iSABR/Raquib-abscopal/853273733430055401-242ac11d-0001-007/'

# patient 106
tumor_file = 'Thr_106_30817-2.igblast.airr.tsv'
nontumor_file = 'Thr_106_30817_1.igblast.airr.tsv'
distal_file = 'Thr_106_30818_1.igblast.airr.tsv'

time_files = [
    '106_0_week.igblast.airr.tsv',
    '106_3_month.igblast.airr.tsv',
    '106_6_month-1.igblast.airr.tsv',
    '106_9_month.igblast.airr.tsv'
]

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

tumor_cdr3 = read_cdr3_dict(dirpath + tumor_file)
nontumor_cdr3 = read_cdr3_dict(dirpath + nontumor_file)
distal_cdr3 = read_cdr3_dict(dirpath + distal_file)
time_cdr3 = []
for f in time_files:
    time_cdr3.append(read_cdr3_dict(dirpath + f))

#
# tumor-associated CDR3s
#

A = set(tumor_cdr3.keys())
B = set(nontumor_cdr3.keys())
venn_cdr3("tumor", A, "non-tumor", B)
TMNT = A - B

i = 0
time_set = []
for f in time_cdr3:
    X = set(f.keys())
    venn_cdr3("TMNT", TMNT, "time " + str(i), X)
    i += 1
    time_set.append(TMNT & X)

X = time_set[0] & time_set[1] & time_set[2] & time_set[3]
venn_cdr3("TMNT", TMNT, "all time", X)
