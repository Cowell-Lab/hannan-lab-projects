import airr
import csv

# AIRR TSV output from igblast
#dirpath = "/work/data/immune/Raquib-iSABR/Raquib-IL2/3401780882221044201-242ac11c-0001-007/"
dirpath = "/work/data/immune/Raquib-iSABR/Raquib-IL2/0dfadecb-9a30-4182-98bd-f07b8437218c-007/"
da_path = '/work/data/immune/Raquib-iSABR/Raquib-IL2/differential_abundance/'

def read_cdr3_dict(filename):
    cdr3 = {}
    data = airr.read_rearrangement(filename)

    for r in data:
        if not r['productive']: continue
        if not r['junction_aa']: continue
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

def check_airr_file(filename, da, cdr3s):
    data = airr.read_rearrangement(filename)
    increase_cnt = 0
    decrease_cnt = 0
    total_cnt = 0
    for row in data:
        for r in da:
            if row['sequence'] == r['sequence']:
                if row['junction_aa'] in cdr3s:
                    print('FOUND IN SET')
                    print(row['sequence'])
                    print(row['junction_aa'], row['productive'], row['duplicate_count'], r['significance'], r['sampleA_abundance'], r['sampleB_abundance'])
                    total_cnt += 1
                    if r['significance'] == 'A > B':
                        decrease_cnt += 1
                    if r['significance'] == 'B > A':
                        increase_cnt += 1
    print('DA FOUND total:', total_cnt)
    print('DA FOUND increase:', increase_cnt)
    print('DA FOUND decrease:', decrease_cnt)

def compare_airr_file(filename1, filename2, da, cdr3s):
    for r in da:
        found1 = False
        data = airr.read_rearrangement(filename1)
        for row in data:
            if row['sequence'] == r['sequence']:
                print(row['sequence'])
                print(row['junction_aa'], row['productive'], row['duplicate_count'], r['significance'], r['sampleA_abundance'], r['sampleB_abundance'])
                if row['junction_aa'] in cdr3s:
                    print('FOUND IN SET')
                found1 = True
        found2 = False
        data = airr.read_rearrangement(filename2)
        for row in data:
            if row['sequence'] == r['sequence']:
                print(row['sequence'])
                print(row['junction_aa'], row['productive'], row['duplicate_count'], r['significance'], r['sampleA_abundance'], r['sampleB_abundance'])
                if row['junction_aa'] in cdr3s:
                    print('FOUND IN SET')
                found2 = True

        if found1 and found2:
            print('FOUND IN BOTH')
        elif found1:
            print('FOUND IN 1')
        elif found2:
            print('FOUND IN 2')


# patient 101
tumor_file_101 = "IL-2_101_PS15-10488_1D-1.igblast.airr.tsv"
nontumor_file_101 = "IL-2_101_PS15-10488_1J-1.igblast.airr.tsv"
pre_file_101 = "101-0wk_TCRB.igblast.airr.tsv"
post_file_101 = "101-9wk_TCRB.igblast.airr.tsv"
da_file_101 = "101-0wk_TCRB_VS_101-9wk_TCRB.differentialAbundance.tsv"

# patient 104
tumor_file_104 = "IL-2_104_12546-13_T-1.igblast.airr.tsv"
nontumor_file_104 = "IL-2_104_12546-3_N-1.igblast.airr.tsv"
pre_file_104 = "104-0wk_TCRB.igblast.airr.tsv"
post_file_104 = "104-9wk_TCRB.igblast.airr.tsv"
da_file_104 = "104-9wk_TCRB_VS_104-0wk_TCRB.differentialAbundance.tsv"

# patient 110
# no tumor for 110, these are bogus files from 104
tumor_file_110 = "IL-2_104_12546-13_T-1.igblast.airr.tsv"
nontumor_file_110 = "IL-2_104_12546-3_N-1.igblast.airr.tsv"
pre_file_110 = "110-0wk_TCRB.igblast.airr.tsv"
post_file_110 = "110-11wk_TCRB.igblast.airr.tsv"
da_file_110 = "110-11wk_TCRB_VS_110-0wk_TCRB.differentialAbundance.tsv"

# patient 115
tumor_file_115 = "IL-2_115_US13-00646-5-T.igblast.airr.tsv"
nontumor_file_115 = "IL-2_115_US13-00646-18-N.igblast.airr.tsv"
pre_file_115 = "JAC-115_0wks.igblast.airr.tsv"
post_file_115 = "JAC-115_11wks.igblast.airr.tsv"
da_file_115 = "JAC-115_11wks_VS_JAC-115_0wks.differentialAbundance.tsv"

# patient 116
tumor_file_116 = "IL-2_116_S13-16236-A6-T.igblast.airr.tsv"
nontumor_file_116 = "IL-2_116_S13-16236-A3-N.igblast.airr.tsv"
pre_file_116 = "116-0wk_TCRB.igblast.airr.tsv"
post_file_116 = "116-11wk_TCRB.igblast.airr.tsv"
da_file_116 = "116-0wk_TCRB_VS_116-11wk_TCRB.differentialAbundance.tsv"

# patient 120
tumor_file_120 = "IL-2_120_25838-5_T.igblast.airr.tsv"
nontumor_file_120 = "IL-2_120_25838-11_N.igblast.airr.tsv"
pre_file_120 = "J-G-120_0wks.igblast.airr.tsv"
post_file_120 = "J-G-120_11wks.igblast.airr.tsv"
da_file_120 = "J-G-120_11wks_VS_J-G-120_0wks.differentialAbundance.tsv"

# patient 122
tumor_file_122 = "IL-2_122_25516-5_T-1.igblast.airr.tsv"
nontumor_file_122 = "IL-2_122_25516-9_N.igblast.airr.tsv"
pre_file_122 = "E-M-122_0wks.igblast.airr.tsv"
post_file_122 = "E-M-122_12wks.igblast.airr.tsv"
da_file_122 = "E-M-122_12wks_VS_E-M-122_0wks.differentialAbundance.tsv"

# patient 128
tumor_file_128 = "IL-2_128_SU16-17418-A3-T.igblast.airr.tsv"
nontumor_file_128 = "IL-2_128_SU16-17418-A1-N.igblast.airr.tsv"
pre_file_128 = "C-H-128_0wks.igblast.airr.tsv"
post_file_128 = "C-H-128_8wks.igblast.airr.tsv"
da_file_128 = "C-H-128_0wks_VS_C-H-128_8wks.differentialAbundance.tsv"

#tumor_file = tumor_file_128
#nontumor_file = nontumor_file_128
#pre_file = pre_file_128
#post_file = post_file_128


#
# tumor-associated CDR3s
#

def new_ta_cdr3(tumor_file, nontumor_file, pre_file, post_file, da_file):
    tumor_cdr3 = read_cdr3_dict(dirpath + tumor_file)
    nontumor_cdr3 = read_cdr3_dict(dirpath + nontumor_file)
    pre_cdr3 = read_cdr3_dict(dirpath + pre_file)
    post_cdr3 = read_cdr3_dict(dirpath + post_file)
    da = read_da_file(da_file)

    A = set(tumor_cdr3.keys())
    B = set(nontumor_cdr3.keys())
    C = set(pre_cdr3.keys())
    D = set(post_cdr3.keys())

    print("DA total: " + str(len(da)))
    print("tumor total: " + str(len(A)))
    print("non-tumor total: " + str(len(B)))
    print("pre-surgery total: " + str(len(C)))
    print("post-surgery total: " + str(len(D)))
    print('')

    X = A & B
    TANT = A & B
    print("TANT shared: " + str(len(X)))
    TMNT = A - B
    print("TMNT unique: " + str(len(TMNT)))
    X = B - A
    print("non-tumor unique: " + str(len(X)))
    print('')

    PMP = D - C
    print("post-surgery PMP unique: " + str(len(PMP)))
    X = TMNT & PMP
    print("TMNT & PMP shared: " + str(len(X)))
    print(X)
    for cdr3 in X:
        print(cdr3 + ',,' + str(tumor_cdr3[cdr3]) + ',' + str(post_cdr3[cdr3]))
    
    check_airr_file(dirpath + post_file, da, X)
    #check_airr_file(dirpath + post_file, da, A)

    return X


def existing_ta_cdr3(tumor_file, nontumor_file, pre_file, post_file, da_file):
    tumor_cdr3 = read_cdr3_dict(dirpath + tumor_file)
    nontumor_cdr3 = read_cdr3_dict(dirpath + nontumor_file)
    pre_cdr3 = read_cdr3_dict(dirpath + pre_file)
    post_cdr3 = read_cdr3_dict(dirpath + post_file)
    da = read_da_file(da_file)

    A = set(tumor_cdr3.keys())
    B = set(nontumor_cdr3.keys())
    C = set(pre_cdr3.keys())
    D = set(post_cdr3.keys())

    print("DA total: " + str(len(da)))
    print("tumor total: " + str(len(A)))
    print("non-tumor total: " + str(len(B)))
    print("pre-surgery total: " + str(len(C)))
    print("post-surgery total: " + str(len(D)))
    print('')

    TMNT = A - B
    PAP = C & D
    print("PAP shared: " + str(len(PAP)))
    X = TMNT & PAP
    print("TMNT & PAP shared: " + str(len(X)))
    print(X)
    for cdr3 in X:
        print(cdr3 + ',,' + str(tumor_cdr3[cdr3]) + ',' + str(pre_cdr3[cdr3]) + ',' + str(post_cdr3[cdr3]))

    check_airr_file(dirpath + pre_file, da, X)
    #check_airr_file(dirpath + post_file, da, X)
    #compare_airr_file(dirpath + pre_file, dirpath + post_file, da, X)

    return X


def old_ta_cdr3(tumor_file, nontumor_file, pre_file, post_file, da_file):
    tumor_cdr3 = read_cdr3_dict(dirpath + tumor_file)
    nontumor_cdr3 = read_cdr3_dict(dirpath + nontumor_file)
    pre_cdr3 = read_cdr3_dict(dirpath + pre_file)
    post_cdr3 = read_cdr3_dict(dirpath + post_file)
    da = read_da_file(da_file)

    A = set(tumor_cdr3.keys())
    B = set(nontumor_cdr3.keys())
    C = set(pre_cdr3.keys())
    D = set(post_cdr3.keys())

    print("DA total: " + str(len(da)))
    print("tumor total: " + str(len(A)))
    print("non-tumor total: " + str(len(B)))
    print("pre-surgery total: " + str(len(C)))
    print("post-surgery total: " + str(len(D)))
    print('')

    TMNT = A - B
    PRE = C - D
    print("pre-surgery unique: " + str(len(PRE)))
    X = TMNT & PRE
    print("TMNT & PRE shared: " + str(len(X)))
    print(X)
    for cdr3 in X:
        print(cdr3 + ',,' + str(tumor_cdr3[cdr3]) + ',' + str(pre_cdr3[cdr3]))

    check_airr_file(dirpath + pre_file, da, X)
    #check_airr_file(dirpath + post_file, da, X)
    #compare_airr_file(dirpath + pre_file, dirpath + post_file, da, X)

    return X

#
# sharing matrix
#
def sharing_matrix(cdr3_sets):
    for i in range(len(cdr3_sets)):
        for j in range(len(cdr3_sets)):
            A = cdr3_sets[i]
            B = cdr3_sets[j]
            C = A & B
            print(i, j, len(C))
            if len(C) > 0:
                print(C)

def gather_new_cdr3_sets():
    new_cdr3_sets = []
    print("*****")
    print("patient 101")
    print("*****")
    new_cdr3_sets.append(new_ta_cdr3(tumor_file_101, nontumor_file_101, pre_file_101, post_file_101, da_file_101))
    print("*****")
    print("patient 104")
    print("*****")
    new_cdr3_sets.append(new_ta_cdr3(tumor_file_104, nontumor_file_104, pre_file_104, post_file_104, da_file_104))
    print("*****")
    print("patient 115")
    print("*****")
    new_cdr3_sets.append(new_ta_cdr3(tumor_file_115, nontumor_file_115, pre_file_115, post_file_115, da_file_115))
    print("*****")
    print("patient 116")
    print("*****")
    new_cdr3_sets.append(new_ta_cdr3(tumor_file_116, nontumor_file_116, pre_file_116, post_file_116, da_file_116))
    print("*****")
    print("patient 120")
    print("*****")
    new_cdr3_sets.append(new_ta_cdr3(tumor_file_120, nontumor_file_120, pre_file_120, post_file_120, da_file_120))
    print("*****")
    print("patient 122")
    print("*****")
    new_cdr3_sets.append(new_ta_cdr3(tumor_file_122, nontumor_file_122, pre_file_122, post_file_122, da_file_122))
    print("*****")
    print("patient 128")
    print("*****")
    new_cdr3_sets.append(new_ta_cdr3(tumor_file_128, nontumor_file_128, pre_file_128, post_file_128, da_file_128))
    return new_cdr3_sets

def gather_existing_cdr3_sets():
    existing_cdr3_sets = []
    print("*****")
    print("patient 101")
    print("*****")
    existing_cdr3_sets.append(existing_ta_cdr3(tumor_file_101, nontumor_file_101, pre_file_101, post_file_101, da_file_101))
    print("*****")
    print("patient 104")
    print("*****")
    existing_cdr3_sets.append(existing_ta_cdr3(tumor_file_104, nontumor_file_104, pre_file_104, post_file_104, da_file_104))
    print("*****")
    print("patient 115")
    print("*****")
    existing_cdr3_sets.append(existing_ta_cdr3(tumor_file_115, nontumor_file_115, pre_file_115, post_file_115, da_file_115))
    print("*****")
    print("patient 116")
    print("*****")
    existing_cdr3_sets.append(existing_ta_cdr3(tumor_file_116, nontumor_file_116, pre_file_116, post_file_116, da_file_116))
    print("*****")
    print("patient 120")
    print("*****")
    existing_cdr3_sets.append(existing_ta_cdr3(tumor_file_120, nontumor_file_120, pre_file_120, post_file_120, da_file_120))
    print("*****")
    print("patient 122")
    print("*****")
    existing_cdr3_sets.append(existing_ta_cdr3(tumor_file_122, nontumor_file_122, pre_file_122, post_file_122, da_file_122))
    print("*****")
    print("patient 128")
    print("*****")
    existing_cdr3_sets.append(existing_ta_cdr3(tumor_file_128, nontumor_file_128, pre_file_128, post_file_128, da_file_128))
    return existing_cdr3_sets

def gather_old_cdr3_sets():
    old_cdr3_sets = []
    print("*****")
    print("patient 101")
    print("*****")
    old_cdr3_sets.append(old_ta_cdr3(tumor_file_101, nontumor_file_101, pre_file_101, post_file_101, da_file_101))
    print("*****")
    print("patient 104")
    print("*****")
    old_cdr3_sets.append(old_ta_cdr3(tumor_file_104, nontumor_file_104, pre_file_104, post_file_104, da_file_104))
    print("*****")
    print("patient 115")
    print("*****")
    old_cdr3_sets.append(old_ta_cdr3(tumor_file_115, nontumor_file_115, pre_file_115, post_file_115, da_file_115))
    print("*****")
    print("patient 116")
    print("*****")
    old_cdr3_sets.append(old_ta_cdr3(tumor_file_116, nontumor_file_116, pre_file_116, post_file_116, da_file_116))
    print("*****")
    print("patient 120")
    print("*****")
    old_cdr3_sets.append(old_ta_cdr3(tumor_file_120, nontumor_file_120, pre_file_120, post_file_120, da_file_120))
    print("*****")
    print("patient 122")
    print("*****")
    old_cdr3_sets.append(old_ta_cdr3(tumor_file_122, nontumor_file_122, pre_file_122, post_file_122, da_file_122))
    print("*****")
    print("patient 128")
    print("*****")
    old_cdr3_sets.append(old_ta_cdr3(tumor_file_128, nontumor_file_128, pre_file_128, post_file_128, da_file_128))
    return old_cdr3_sets


#####
#
# analysis
#

existing_ta_cdr3(tumor_file_110, nontumor_file_110, pre_file_110, post_file_110, da_file_110)

#new_cdr3_sets = gather_new_cdr3_sets()
#sharing_matrix(new_cdr3_sets)

#existing_cdr3_sets = gather_existing_cdr3_sets()
#sharing_matrix(existing_cdr3_sets)

#old_cdr3_sets = gather_old_cdr3_sets()
#sharing_matrix(old_cdr3_sets)
