import airr
import operator

igblast_dir = '/work/data/immune/Raquib-iSABR/Raquib-abscopal/853273733430055401-242ac11d-0001-007/'

# patient 101
#tumor_file = 'Thr_101_25516-1.igblast.airr.tsv'
#nontumor_file = 'Thr_101_25516_9.igblast.airr.tsv'
#distal_file = 'Thr_101_25519-1.igblast.airr.tsv'

# patient 103
#tumor_file = 'Thr_103_26277_1.igblast.airr.tsv'
#nontumor_file = 'Thr_103_26277_3.igblast.airr.tsv'
#distal_file = 'Thr_103_26276-1.igblast.airr.tsv'

# patient 104
tumor_file = '24673-6.igblast.airr.tsv'
nontumor_file = '24563-5.igblast.airr.tsv'
distal_file = '26473-1.igblast.airr.tsv'

# patient 105
#tumor_file = 'Thr_105_28473_2H.igblast.airr.tsv'
#nontumor_file = 'Thr_105_28473_3.igblast.airr.tsv'
#distal_file = 'Thr_105_28474-4.igblast.airr.tsv'

# patient 106
#tumor_file = 'Thr_106_30817-2.igblast.airr.tsv'
#nontumor_file = 'Thr_106_30817_1.igblast.airr.tsv'
#distal_file = 'Thr_106_30818_1.igblast.airr.tsv'

# patient 107
#tumor_file = 'Thr_107_31626-3.igblast.airr.tsv'
#nontumor_file = 'Thr_107_31626_1.igblast.airr.tsv'
#distal_file = 'Thr_107_31627-1.igblast.airr.tsv'

tumor_cdr3s = {}
nontumor_cdr3s = {}
distal_cdr3s = {}

reader = airr.read_rearrangement(igblast_dir + tumor_file)
for row in reader:
    if not row['productive']: continue
    if tumor_cdr3s.get(row['junction_aa']) is None:
        tumor_cdr3s[row['junction_aa']] = row['duplicate_count']
    else:
        tumor_cdr3s[row['junction_aa']] += row['duplicate_count']

reader = airr.read_rearrangement(igblast_dir + nontumor_file)
for row in reader:
    if not row['productive']: continue
    if nontumor_cdr3s.get(row['junction_aa']) is None:
        nontumor_cdr3s[row['junction_aa']] = row['duplicate_count']
    else:
        nontumor_cdr3s[row['junction_aa']] += row['duplicate_count']

reader = airr.read_rearrangement(igblast_dir + distal_file)
for row in reader:
    if not row['productive']: continue
    if distal_cdr3s.get(row['junction_aa']) is None:
        distal_cdr3s[row['junction_aa']] = row['duplicate_count']
    else:
        distal_cdr3s[row['junction_aa']] += row['duplicate_count']

print('tumor total',len(tumor_cdr3s))
print('nontumor total',len(nontumor_cdr3s))
print('distal total',len(distal_cdr3s))

tumor_set = set(tumor_cdr3s.keys())
nontumor_set = set(nontumor_cdr3s.keys())
distal_set = set(distal_cdr3s.keys())

all_share = tumor_set & nontumor_set & distal_set
print('all intersect', len(all_share))
X = tumor_set & nontumor_set - all_share
print('tumor,nontumor intersect',len(X))
X = tumor_set & distal_set - all_share
print('tumor,distal intersect',len(X))
X = nontumor_set & distal_set - all_share
print('nontumor,distal intersect',len(X))
X = tumor_set - nontumor_set - distal_set
print('tumor only',len(X))
X = nontumor_set - tumor_set - distal_set
print('nontumor only',len(X))
X = distal_set - nontumor_set - tumor_set
print('distal only',len(X))

sorted_x = sorted(tumor_cdr3s.items(), key=operator.itemgetter(1), reverse=True)
print('top tumor',sorted_x[0:9])
sorted_x = sorted(nontumor_cdr3s.items(), key=operator.itemgetter(1), reverse=True)
print('top nontumor',sorted_x[0:9])
sorted_x = sorted(distal_cdr3s.items(), key=operator.itemgetter(1), reverse=True)
print('top distal',sorted_x[0:9])

