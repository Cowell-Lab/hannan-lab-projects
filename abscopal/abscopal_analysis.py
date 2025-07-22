import airr

#dirpath = "105060807991225880-242ac11b-0001-007/"
dirpath = "2403978270247808535-242ac11b-0001-007/"

tumor_file = "24673-6.igblast.airr.tsv"
nontumor_file = "24563-5.igblast.airr.tsv"
lung_file = "022015-058-7.igblast.airr.tsv"
time_files = [
    "O_weeks.igblast.airr.tsv",
    "5_weeks.igblast.airr.tsv",
    "5_weeks_2.igblast.airr.tsv",
    "8_weeks.igblast.airr.tsv",
    "17_weeks.igblast.airr.tsv",
    "7_months.igblast.airr.tsv",
    "30_months.igblast.airr.tsv"
    ]

def read_cdr3_dict(filename):
    cdr3 = {}
    data = airr.read_rearrangement(open(filename, 'r'))

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
lung_cdr3 = read_cdr3_dict(dirpath + lung_file)
time_cdr3 = []
for f in time_files:
    time_cdr3.append(read_cdr3_dict(dirpath + f))

#
# tumor-associated CDR3s
#

A = set(tumor_cdr3.keys())
B = set(nontumor_cdr3.keys())
C = set(lung_cdr3.keys())

print("tumor total: " + str(len(A)))
print("non-tumor total: " + str(len(B)))
print("lung biopsy total: " + str(len(C)))
X = A & B
TANT = A & B
print("shared: " + str(len(X)))
TMNT = A - B
print("tumor unique: " + str(len(TMNT)))
X = B - A
print("non-tumor unique: " + str(len(X)))
X = TMNT & C
print("shared: " + str(len(X)))
X = TMNT - C
print("tumor unique: " + str(len(X)))
X = C - TMNT
print("lung unique: " + str(len(X)))

venn_cdr3("tumor", A, "non-tumor", B)
venn_cdr3("TMNT", TMNT, "lung", C)
venn_cdr3("non-tumor", B, "lung", C)

i = 0
time_set = []
for f in time_cdr3:
    X = set(f.keys())
    venn_cdr3("TMNT", TMNT, "time " + str(i), X)
    i += 1
    time_set.append(TMNT & X)

X = time_set[0] & time_set[1] & time_set[2] & time_set[3] & time_set[4] & time_set[5] & time_set[6]
venn_cdr3("TMNT", TMNT, "all time", X)

X = TMNT & X
for f in X:
    print(f +','+ str(time_cdr3[0][f]) +','+ str(time_cdr3[1][f])
          +','+ str(time_cdr3[2][f]) +','+ str(time_cdr3[3][f]) +','+ str(time_cdr3[4][f]) +','+ str(time_cdr3[5][f])
          +','+ str(time_cdr3[6][f]) +','+ str(tumor_cdr3[f]))

venn_cdr3("TMNT, all time", X, "lung", C)
X = X & C
for f in X: print(f +','+ str(tumor_cdr3[f]) +','+ str(lung_cdr3[f]))

i = 0
time_set = []
for f in time_cdr3:
    X = set(f.keys())
    venn_cdr3("T&NT", TANT, "time " + str(i), X)
    i += 1
    time_set.append(TANT & X)

X = time_set[0] & time_set[1] & time_set[2] & time_set[3] & time_set[4] & time_set[5] & time_set[6]
venn_cdr3("T&NT", TANT, "all time", X)
venn_cdr3("T&NT", TANT, "lung", C)

X = TANT & X
for f in X:
    print(f +','+ str(time_cdr3[0][f]) +','+ str(time_cdr3[1][f])
          +','+ str(time_cdr3[2][f]) +','+ str(time_cdr3[3][f]) +','+ str(time_cdr3[4][f]) +','+ str(time_cdr3[5][f])
          +','+ str(time_cdr3[6][f])
          +','+ str(tumor_cdr3[f]) +','+ str(nontumor_cdr3[f]))

venn_cdr3("TANT, all time", X, "lung", C)
X = X & C
for f in X: print(f +','+ str(tumor_cdr3[f]) +','+ str(lung_cdr3[f]))

time_set = []
for f in time_cdr3:
    time_set.append(set(f.keys()))

T5 = time_set[1] | time_set[2]
venn_cdr3("0 weeks", time_set[0], "5 weeks", T5)

X = T5 - time_set[0]
venn_cdr3("5 weeks", X, "TMNT", TMNT)
venn_cdr3("5 weeks", X, "T&NT", TANT)

X = X & TMNT
for f in X:
    n1 = ''
    n2 = ''
    if time_cdr3[1].get(f): n1 = str(time_cdr3[1][f])
    if time_cdr3[2].get(f): n2 = str(time_cdr3[2][f])
    print(f +','+ n1 +','+ n2 +','+ str(tumor_cdr3[f]))

T1P = time_set[1] | time_set[2] | time_set[3] | time_set[4] | time_set[5] | time_set[6]
venn_cdr3("0 weeks", time_set[0], "after 0", T1P)

T0P = time_set[0] | time_set[1] | time_set[2] | time_set[3] | time_set[4] | time_set[5] | time_set[6]
venn_cdr3("any time", T0P, "TMNT", TMNT)
venn_cdr3("any time", T0P, "T&NT", TANT)

T1P = T0P - time_set[0]
venn_cdr3("only after 0", T1P, "TMNT", TMNT)
venn_cdr3("only after 0", T1P, "T&NT", TANT)

T2P = (T1P - time_set[1]) - time_set[2]
venn_cdr3("only after 5", T2P, "TMNT", TMNT)
venn_cdr3("only after 5", T2P, "T&NT", TANT)

T3P = T2P - time_set[3]
venn_cdr3("only after 8", T3P, "TMNT", TMNT)
venn_cdr3("only after 8", T3P, "T&NT", TANT)

T4P = T3P - time_set[4]
venn_cdr3("only after 17", T4P, "TMNT", TMNT)
venn_cdr3("only after 17", T4P, "T&NT", TANT)

T5P = T4P - time_set[5]
venn_cdr3("only 7 months", T5P, "TMNT", TMNT)
venn_cdr3("only 7 months", T5P, "T&NT", TANT)

X = T5P & TMNT
for f in X:
    n0 = ''
    n1 = ''
    n2 = ''
    n3 = ''
    n4 = ''
    n5 = ''
    n6 = ''
    if time_cdr3[0].get(f): n0 = str(time_cdr3[0][f])
    if time_cdr3[1].get(f): n1 = str(time_cdr3[1][f])
    if time_cdr3[2].get(f): n2 = str(time_cdr3[2][f])
    if time_cdr3[3].get(f): n3 = str(time_cdr3[3][f])
    if time_cdr3[4].get(f): n4 = str(time_cdr3[4][f])
    if time_cdr3[5].get(f): n5 = str(time_cdr3[5][f])
    if time_cdr3[6].get(f): n6 = str(time_cdr3[6][f])
    print(f +','+ n0 +','+ n1 +','+ n2 +','+ n3 +','+ n4 +','+ n5 +','+ n6 +',,'+ str(tumor_cdr3[f]))

print('')
X = (T4P - T5P) & TMNT
for f in X:
    n0 = ''
    n1 = ''
    n2 = ''
    n3 = ''
    n4 = ''
    n5 = ''
    n6 = ''
    if time_cdr3[0].get(f): n0 = str(time_cdr3[0][f])
    if time_cdr3[1].get(f): n1 = str(time_cdr3[1][f])
    if time_cdr3[2].get(f): n2 = str(time_cdr3[2][f])
    if time_cdr3[3].get(f): n3 = str(time_cdr3[3][f])
    if time_cdr3[4].get(f): n4 = str(time_cdr3[4][f])
    if time_cdr3[5].get(f): n5 = str(time_cdr3[5][f])
    if time_cdr3[6].get(f): n6 = str(time_cdr3[6][f])
    print(f +','+ n0 +','+ n1 +','+ n2 +','+ n3 +','+ n4 +','+ n5 +','+ n6 +',,'+ str(tumor_cdr3[f]))

print('')
X = (T3P - T4P - T5P) & TMNT
for f in X:
    n0 = ''
    n1 = ''
    n2 = ''
    n3 = ''
    n4 = ''
    n5 = ''
    n6 = ''
    if time_cdr3[0].get(f): n0 = str(time_cdr3[0][f])
    if time_cdr3[1].get(f): n1 = str(time_cdr3[1][f])
    if time_cdr3[2].get(f): n2 = str(time_cdr3[2][f])
    if time_cdr3[3].get(f): n3 = str(time_cdr3[3][f])
    if time_cdr3[4].get(f): n4 = str(time_cdr3[4][f])
    if time_cdr3[5].get(f): n5 = str(time_cdr3[5][f])
    if time_cdr3[6].get(f): n6 = str(time_cdr3[6][f])
    print(f +','+ n0 +','+ n1 +','+ n2 +','+ n3 +','+ n4 +','+ n5 +','+ n6 +',,'+ str(tumor_cdr3[f]))

print('')
X = (T2P - T3P - T4P - T5P) & TMNT
for f in X:
    n0 = ''
    n1 = ''
    n2 = ''
    n3 = ''
    n4 = ''
    n5 = ''
    n6 = ''
    if time_cdr3[0].get(f): n0 = str(time_cdr3[0][f])
    if time_cdr3[1].get(f): n1 = str(time_cdr3[1][f])
    if time_cdr3[2].get(f): n2 = str(time_cdr3[2][f])
    if time_cdr3[3].get(f): n3 = str(time_cdr3[3][f])
    if time_cdr3[4].get(f): n4 = str(time_cdr3[4][f])
    if time_cdr3[5].get(f): n5 = str(time_cdr3[5][f])
    if time_cdr3[6].get(f): n6 = str(time_cdr3[6][f])
    print(f +','+ n0 +','+ n1 +','+ n2 +','+ n3 +','+ n4 +','+ n5 +','+ n6 +',,'+ str(tumor_cdr3[f]))

print('')
X = (T1P - T2P - T3P - T4P - T5P) & TMNT
for f in X:
    n0 = ''
    n1 = ''
    n2 = ''
    n3 = ''
    n4 = ''
    n5 = ''
    n6 = ''
    if time_cdr3[0].get(f): n0 = str(time_cdr3[0][f])
    if time_cdr3[1].get(f): n1 = str(time_cdr3[1][f])
    if time_cdr3[2].get(f): n2 = str(time_cdr3[2][f])
    if time_cdr3[3].get(f): n3 = str(time_cdr3[3][f])
    if time_cdr3[4].get(f): n4 = str(time_cdr3[4][f])
    if time_cdr3[5].get(f): n5 = str(time_cdr3[5][f])
    if time_cdr3[6].get(f): n6 = str(time_cdr3[6][f])
    print(f +','+ n0 +','+ n1 +','+ n2 +','+ n3 +','+ n4 +','+ n5 +','+ n6 +',,'+ str(tumor_cdr3[f]))

print('')
Y = time_set[0] & time_set[1] & time_set[2] & time_set[3] & time_set[4] & time_set[5] & time_set[6]
Y = Y & TMNT
X = (T0P - T1P - T2P - T3P - T4P - T5P - Y) & TMNT
for f in X:
    n0 = ''
    n1 = ''
    n2 = ''
    n3 = ''
    n4 = ''
    n5 = ''
    n6 = ''
    if time_cdr3[0].get(f): n0 = str(time_cdr3[0][f])
    if time_cdr3[1].get(f): n1 = str(time_cdr3[1][f])
    if time_cdr3[2].get(f): n2 = str(time_cdr3[2][f])
    if time_cdr3[3].get(f): n3 = str(time_cdr3[3][f])
    if time_cdr3[4].get(f): n4 = str(time_cdr3[4][f])
    if time_cdr3[5].get(f): n5 = str(time_cdr3[5][f])
    if time_cdr3[6].get(f): n6 = str(time_cdr3[6][f])
    print(f +','+ n0 +','+ n1 +','+ n2 +','+ n3 +','+ n4 +','+ n5 +','+ n6 +',,'+ str(tumor_cdr3[f]))

print('****')

X = T5P & TANT
for f in X:
    n0 = ''
    n1 = ''
    n2 = ''
    n3 = ''
    n4 = ''
    n5 = ''
    n6 = ''
    if time_cdr3[0].get(f): n0 = str(time_cdr3[0][f])
    if time_cdr3[1].get(f): n1 = str(time_cdr3[1][f])
    if time_cdr3[2].get(f): n2 = str(time_cdr3[2][f])
    if time_cdr3[3].get(f): n3 = str(time_cdr3[3][f])
    if time_cdr3[4].get(f): n4 = str(time_cdr3[4][f])
    if time_cdr3[5].get(f): n5 = str(time_cdr3[5][f])
    if time_cdr3[6].get(f): n6 = str(time_cdr3[6][f])
    print(f +','+ n0 +','+ n1 +','+ n2 +','+ n3 +','+ n4 +','+ n5 +','+ n6 +',,'+ str(tumor_cdr3[f]) +','+ str(nontumor_cdr3[f]))

print('')
X = (T4P - T5P) & TANT
for f in X:
    n0 = ''
    n1 = ''
    n2 = ''
    n3 = ''
    n4 = ''
    n5 = ''
    n6 = ''
    if time_cdr3[0].get(f): n0 = str(time_cdr3[0][f])
    if time_cdr3[1].get(f): n1 = str(time_cdr3[1][f])
    if time_cdr3[2].get(f): n2 = str(time_cdr3[2][f])
    if time_cdr3[3].get(f): n3 = str(time_cdr3[3][f])
    if time_cdr3[4].get(f): n4 = str(time_cdr3[4][f])
    if time_cdr3[5].get(f): n5 = str(time_cdr3[5][f])
    if time_cdr3[6].get(f): n6 = str(time_cdr3[6][f])
    print(f +','+ n0 +','+ n1 +','+ n2 +','+ n3 +','+ n4 +','+ n5 +','+ n6 +',,'+ str(tumor_cdr3[f]) +','+ str(nontumor_cdr3[f]))

print('')
X = (T3P - T4P - T5P) & TANT
for f in X:
    n0 = ''
    n1 = ''
    n2 = ''
    n3 = ''
    n4 = ''
    n5 = ''
    n6 = ''
    if time_cdr3[0].get(f): n0 = str(time_cdr3[0][f])
    if time_cdr3[1].get(f): n1 = str(time_cdr3[1][f])
    if time_cdr3[2].get(f): n2 = str(time_cdr3[2][f])
    if time_cdr3[3].get(f): n3 = str(time_cdr3[3][f])
    if time_cdr3[4].get(f): n4 = str(time_cdr3[4][f])
    if time_cdr3[5].get(f): n5 = str(time_cdr3[5][f])
    if time_cdr3[6].get(f): n6 = str(time_cdr3[6][f])
    print(f +','+ n0 +','+ n1 +','+ n2 +','+ n3 +','+ n4 +','+ n5 +','+ n6 +',,'+ str(tumor_cdr3[f]) +','+ str(nontumor_cdr3[f]))

print('')
X = (T2P - T3P - T4P - T5P) & TANT
for f in X:
    n0 = ''
    n1 = ''
    n2 = ''
    n3 = ''
    n4 = ''
    n5 = ''
    n6 = ''
    if time_cdr3[0].get(f): n0 = str(time_cdr3[0][f])
    if time_cdr3[1].get(f): n1 = str(time_cdr3[1][f])
    if time_cdr3[2].get(f): n2 = str(time_cdr3[2][f])
    if time_cdr3[3].get(f): n3 = str(time_cdr3[3][f])
    if time_cdr3[4].get(f): n4 = str(time_cdr3[4][f])
    if time_cdr3[5].get(f): n5 = str(time_cdr3[5][f])
    if time_cdr3[6].get(f): n6 = str(time_cdr3[6][f])
    print(f +','+ n0 +','+ n1 +','+ n2 +','+ n3 +','+ n4 +','+ n5 +','+ n6 +',,'+ str(tumor_cdr3[f]) +','+ str(nontumor_cdr3[f]))

print('')
X = (T1P - T2P - T3P - T4P - T5P) & TANT
for f in X:
    n0 = ''
    n1 = ''
    n2 = ''
    n3 = ''
    n4 = ''
    n5 = ''
    n6 = ''
    if time_cdr3[0].get(f): n0 = str(time_cdr3[0][f])
    if time_cdr3[1].get(f): n1 = str(time_cdr3[1][f])
    if time_cdr3[2].get(f): n2 = str(time_cdr3[2][f])
    if time_cdr3[3].get(f): n3 = str(time_cdr3[3][f])
    if time_cdr3[4].get(f): n4 = str(time_cdr3[4][f])
    if time_cdr3[5].get(f): n5 = str(time_cdr3[5][f])
    if time_cdr3[6].get(f): n6 = str(time_cdr3[6][f])
    print(f +','+ n0 +','+ n1 +','+ n2 +','+ n3 +','+ n4 +','+ n5 +','+ n6 +',,'+ str(tumor_cdr3[f]) +','+ str(nontumor_cdr3[f]))

print('')
Y = time_set[0] & time_set[1] & time_set[2] & time_set[3] & time_set[4] & time_set[5] & time_set[6]
Y = Y & TANT
X = (T0P - T1P - T2P - T3P - T4P - T5P - Y) & TANT
for f in X:
    n0 = ''
    n1 = ''
    n2 = ''
    n3 = ''
    n4 = ''
    n5 = ''
    n6 = ''
    if time_cdr3[0].get(f): n0 = str(time_cdr3[0][f])
    if time_cdr3[1].get(f): n1 = str(time_cdr3[1][f])
    if time_cdr3[2].get(f): n2 = str(time_cdr3[2][f])
    if time_cdr3[3].get(f): n3 = str(time_cdr3[3][f])
    if time_cdr3[4].get(f): n4 = str(time_cdr3[4][f])
    if time_cdr3[5].get(f): n5 = str(time_cdr3[5][f])
    if time_cdr3[6].get(f): n6 = str(time_cdr3[6][f])
    print(f +','+ n0 +','+ n1 +','+ n2 +','+ n3 +','+ n4 +','+ n5 +','+ n6 +',,'+ str(tumor_cdr3[f]) +','+ str(nontumor_cdr3[f]))

#
# radiation-associated CDR3s
#

# 46556 * 0.01% ~ 5
# 46556 * 0.05% ~ 23
baseline_threshold = 5
baseline_above = {}
baseline_below = {}
for cdr3 in time_cdr3[0]:
    if time_cdr3[0][cdr3] < baseline_threshold:
        baseline_below[cdr3] = time_cdr3[0][cdr3]
    else:
        baseline_above[cdr3] = time_cdr3[0][cdr3]

# 14713 * 0.05% ~ 7
# 14713 * 0.01% ~ 2
s1_threshold = 7
#s1_threshold = 2
s1_above = {}
s1_below = {}
for cdr3 in time_cdr3[1]:
    if time_cdr3[1][cdr3] < s1_threshold:
        s1_below[cdr3] = time_cdr3[1][cdr3]
    else:
        s1_above[cdr3] = time_cdr3[1][cdr3]
# 25978 * 0.05% ~ 13
# 25978 * 0.01% ~ 3
s2_threshold = 13
s2_above = {}
s2_below = {}
for cdr3 in time_cdr3[2]:
    if time_cdr3[2][cdr3] < s2_threshold:
        s2_below[cdr3] = time_cdr3[2][cdr3]
    else:
        s2_above[cdr3] = time_cdr3[2][cdr3]

A = set(s1_above.keys()) | set(s2_above.keys())
B = set(baseline_above.keys())
C = A - B
D1 = {}
D2 = {}
for cdr3 in C:
    if s1_above.get(cdr3):
        D1[cdr3] = s1_above[cdr3]
for cdr3 in C:
    if s2_above.get(cdr3):
        D2[cdr3] = s2_above[cdr3]

print('****')
print('radiation-associated CDR3s')
print('****')

for f in C:
    n0 = ''
    n1 = ''
    n2 = ''
    n3 = ''
    n4 = ''
    n5 = ''
    n6 = ''
    n7 = ''
    n8 = ''
    n9 = ''
    if time_cdr3[0].get(f): n0 = str(time_cdr3[0][f])
    if time_cdr3[1].get(f): n1 = str(time_cdr3[1][f])
    if time_cdr3[2].get(f): n2 = str(time_cdr3[2][f])
    if time_cdr3[3].get(f): n3 = str(time_cdr3[3][f])
    if time_cdr3[4].get(f): n4 = str(time_cdr3[4][f])
    if time_cdr3[5].get(f): n5 = str(time_cdr3[5][f])
    if time_cdr3[6].get(f): n6 = str(time_cdr3[6][f])
    if tumor_cdr3.get(f): n7 = str(tumor_cdr3[f])
    if nontumor_cdr3.get(f): n8 = str(nontumor_cdr3[f])
    if lung_cdr3.get(f): n9 = str(lung_cdr3[f])
    print(f +','+ n0 +','+ n1 +','+ n2 +','+ n3 +','+ n4 +','+ n5 +','+ n6 +',,'+ n7 +','+ n8 +','+ n9)
