
def within_group():
    #fin = open('/Users/s166813/Projects/data/immune/Raquib-iSABR/Raquib-IL2/3069485717384326680-242ac11b-0001-007/cdr3_sharing_data/sampleGroup5_cdr3_aa_sharing.tsv', 'r')
    fin = open('/Users/s166813/Projects/data/immune/Raquib-iSABR/Raquib-IL2/3069485717384326680-242ac11b-0001-007/cdr3_sharing_data/sampleGroup1_cdr3_aa_sharing.tsv', 'r')

    p1cnt = {}
    p1cdr3 = {}

    first = True
    while True:
        line = fin.readline()
        if first:
            first = False
            continue
        if not line:
            break

        fields = line.split('\t')
        if p1cnt.get(fields[1]) is None:
            p1cnt[fields[1]] = 1
            p1cdr3[fields[1]] = [ fields[0] ]
        else:
            p1cnt[fields[1]] += 1
            p1cdr3[fields[1]].append(fields[0])

        if fields[1] == '5690109703415534056-242ac119-0001-012,5623709509019374056-242ac119-0001-012,6106936279492334056-242ac119-0001-012':
            print(fields)
        if fields[1] == '5658455794444014056-242ac119-0001-012,6076270212998894056-242ac119-0001-012,5583852212512494056-242ac119-0001-012':
            print(fields)
        if fields[1] == '5658455794444014056-242ac119-0001-012,5583852212512494056-242ac119-0001-012':
            print(fields)
#        if fields[1] == '6076270212998894056-242ac119-0001-012,5583852212512494056-242ac119-0001-012':
#            print(fields)
#        if fields[1] == '5658455794444014056-242ac119-0001-012,6076270212998894056-242ac119-0001-012':
#            print(fields)
    print(p1cnt)
    #print(p1cdr3['5690109703415534056-242ac119-0001-012,5623709509019374056-242ac119-0001-012,6106936279492334056-242ac119-0001-012'])
    print(p1cdr3['5658455794444014056-242ac119-0001-012,6076270212998894056-242ac119-0001-012,5583852212512494056-242ac119-0001-012'])


def between_group():
    fin = open('/Users/s166813/Projects/data/immune/Raquib-iSABR/Raquib-IL2/3069485717384326680-242ac11b-0001-007/cdr3_sharing_data/group_sampleGroup_comparison_cdr3_aa_sharing.tsv', 'r')

    found = False
    while True:
        line = fin.readline()
        if not line:
            break

        line = line.replace('\n','')
        fields = line.split('\t')
        if fields[0] == 'GROUPS':
            if fields[1] == 'sampleGroup1' and fields[2] == 'sampleGroup5':
                print(fields)
                found = True
                continue
            else:
                found = False
        else:
            if not found:
                continue
            if fields[0] == 'CDR3 AA (imgt)':
                continue
            if len(fields) < 7:
                continue
            if int(fields[2]) == 3 and int(fields[6]) == 3:
                print(fields)


def diff_between_group():
    fin = open('/Users/s166813/Projects/data/immune/Raquib-iSABR/Raquib-IL2/3069485717384326680-242ac11b-0001-007/cdr3_sharing_data/group_sampleGroup_diff_cdr3_aa_sharing.tsv', 'r')

    found = False
    while True:
        line = fin.readline()
        if not line:
            break

        line = line.replace('\n','')
        fields = line.split('\t')
        if fields[0] == 'GROUPS':
            if fields[1] == 'sampleGroup1' and fields[2] == 'sampleGroup5':
                print(fields)
                found = True
                continue
            if fields[1] == 'sampleGroup5' and fields[2] == 'sampleGroup1':
                print(fields)
                found = True
                continue
            found = False
        else:
            if not found:
                continue
            if fields[0] == 'CDR3 AA (imgt)':
                continue
            if int(fields[2]) == 3:
                print(fields)





within_group()
#between_group()
#diff_between_group()
