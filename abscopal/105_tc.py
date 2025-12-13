import airr

dirpath = '~/Projects/HannanLab/abscopal/vdjserver/853273733430055401-242ac11d-0001-007/'

# patient 105
tumor_file = 'Thr_105_28473_2H.igblast.airr.tsv'
nontumor_file = 'Thr_105_28473_3.igblast.airr.tsv'
distal_file = 'Thr_105_28474-4.igblast.airr.tsv'

time_files = [
    '105_0_week.igblast.airr.tsv',
    '105_2_week.igblast.airr.tsv',
    '105_3_month.igblast.airr.tsv',
    '105_6_month.igblast.airr.tsv'
]

# def read_cdr3_dict(filename):
#     cdr3 = {}
#     data = airr.read_rearrangement(filename)

#     for r in data:
#         if not r['productive']: continue
#         if cdr3.get(r['junction_aa']):
#             cdr3[r['junction_aa']] = cdr3[r['junction_aa']] + r['duplicate_count']
#         else:
#             cdr3[r['junction_aa']] = r['duplicate_count']

#     return cdr3

def read_cdr3_dict(filename):
    print(filename)
    cdr3 = {}
    data = pd.read_csv(filename, sep = '\t', usecols = ['productive', 'junction_aa', 'duplicate_count'])
    data = data[data.productive == 'T'].drop("productive", axis=1)
    data = data.groupby(data.junction_aa).aggregate({"duplicate_count": 'sum'}).reset_index()
    data.set_index('junction_aa', inplace = True)
    result_dict = data['duplicate_count'].to_dict()
    return result_dict

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

X = time_set[0] & time_set[2] & time_set[3]
venn_cdr3("TMNT", TMNT, "all time", X)
