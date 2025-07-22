#
# Conversion script for Hannan Lab IL-2 study
#
# First run the standard convert_repertoires.py
# Then run this to do specific cleanup
#

import airr

rep_file = './repertoires.airr.json'
out_file = './IL2.airr.json'

data = airr.load_repertoire(rep_file)
reps = data['Repertoire']

for r in reps:
    # hand coded stuff

    for entry in r['sample']:
        if entry['tissue'] == 'tumor':
            entry['tissue'] = { 'id': 'UBERON:0002113', 'label': 'kidney' }
            entry['sample_type'] = 'Formalin-fixed paraffin-embedded'
            entry['cell_storage'] = False
        elif entry['tissue'] == 'nontumor':
            entry['tissue'] = { 'id': 'UBERON:0002113', 'label': 'kidney' }
            entry['sample_type'] = 'Formalin-fixed paraffin-embedded'
            entry['cell_storage'] = False
        elif entry['tissue'] == 'PBMC':
            entry['tissue'] = { 'id': 'UBERON:0013756', 'label': 'venous blood' }
            entry['cell_storage'] = True
        else:
            print('unknown tissue', entry['tissue'])
        entry['template_class'] = 'DNA'
        entry['library_generation_method'] = 'PCR'
        entry['library_generation_protocol'] = 'Adaptive Biotechnologies'
        entry['library_generation_kit_version'] = 'v3'
        entry['pcr_target'][0]['pcr_target_locus'] = 'TRB'
        entry['cell_subset'] = None
        entry['complete_sequences'] = 'partial'
        entry['physical_linkage'] = 'none'
        entry['sequencing_files']['file_type'] = 'fasta'
        entry['sequencing_files']['read_direction'] = 'forward'
        if '0wk' in entry['sample_id']:
            entry['collection_time_point_relative'] = 0
            entry['collection_time_point_relative_unit'] = { 'id': 'UO:0000034', 'label': 'week' }
        elif '8wk' in entry['sample_id']:
            entry['collection_time_point_relative'] = 8
            entry['collection_time_point_relative_unit'] = { 'id': 'UO:0000034', 'label': 'week' }
        elif '9wk' in entry['sample_id']:
            entry['collection_time_point_relative'] = 9
            entry['collection_time_point_relative_unit'] = { 'id': 'UO:0000034', 'label': 'week' }
        elif '11wk' in entry['sample_id']:
            entry['collection_time_point_relative'] = 11
            entry['collection_time_point_relative_unit'] = { 'id': 'UO:0000034', 'label': 'week' }
        elif '12wk' in entry['sample_id']:
            entry['collection_time_point_relative'] = 12
            entry['collection_time_point_relative_unit'] = { 'id': 'UO:0000034', 'label': 'week' }
        else:
            entry['collection_time_point_relative'] = 0
            entry['collection_time_point_relative_unit'] = { 'id': 'UO:0000034', 'label': 'week' }

    if r['subject']['subject_id'] == '101':
        r['subject']['sex'] = 'male'
        r['subject']['race'] = 'White'
        r['subject']['age_min'] = 43
        r['subject']['age_max'] = 43
        r['subject']['age_unit']['id'] = 'UO:0000036'
        r['subject']['age_unit']['label'] = 'year'
        r['subject']['diagnosis'][0]['study_group_description'] = 'CR'
        r['subject']['diagnosis'][0]['disease_diagnosis']['id'] = 'DOID:4450'
        r['subject']['diagnosis'][0]['disease_diagnosis']['label'] = 'renal cell carcinoma'
    elif r['subject']['subject_id'] == '104':
        r['subject']['sex'] = 'male'
        r['subject']['race'] = 'White'
        r['subject']['age_min'] = 54
        r['subject']['age_max'] = 54
        r['subject']['age_unit']['id'] = 'UO:0000036'
        r['subject']['age_unit']['label'] = 'year'
        r['subject']['diagnosis'][0]['study_group_description'] = 'SD'
        r['subject']['diagnosis'][0]['disease_diagnosis']['id'] = 'DOID:4450'
        r['subject']['diagnosis'][0]['disease_diagnosis']['label'] = 'renal cell carcinoma'
    elif r['subject']['subject_id'] == '110':
        r['subject']['sex'] = 'male'
        r['subject']['race'] = 'White'
        r['subject']['age_min'] = 50
        r['subject']['age_max'] = 50
        r['subject']['age_unit']['id'] = 'UO:0000036'
        r['subject']['age_unit']['label'] = 'year'
        r['subject']['diagnosis'][0]['study_group_description'] = 'PD'
        r['subject']['diagnosis'][0]['disease_diagnosis']['id'] = 'DOID:4450'
        r['subject']['diagnosis'][0]['disease_diagnosis']['label'] = 'renal cell carcinoma'
    elif r['subject']['subject_id'] == '115':
        r['subject']['sex'] = 'female'
        r['subject']['race'] = 'White'
        r['subject']['age_min'] = 48
        r['subject']['age_max'] = 48
        r['subject']['age_unit']['id'] = 'UO:0000036'
        r['subject']['age_unit']['label'] = 'year'
        r['subject']['diagnosis'][0]['study_group_description'] = 'PD'
        r['subject']['diagnosis'][0]['disease_diagnosis']['id'] = 'DOID:4450'
        r['subject']['diagnosis'][0]['disease_diagnosis']['label'] = 'renal cell carcinoma'
    elif r['subject']['subject_id'] == '116':
        r['subject']['sex'] = 'male'
        r['subject']['race'] = 'White'
        r['subject']['age_min'] = 44
        r['subject']['age_max'] = 44
        r['subject']['age_unit']['id'] = 'UO:0000036'
        r['subject']['age_unit']['label'] = 'year'
        r['subject']['diagnosis'][0]['study_group_description'] = 'PD'
        r['subject']['diagnosis'][0]['disease_diagnosis']['id'] = 'DOID:4450'
        r['subject']['diagnosis'][0]['disease_diagnosis']['label'] = 'renal cell carcinoma'
    elif r['subject']['subject_id'] == '120':
        r['subject']['sex'] = 'male'
        r['subject']['race'] = 'White'
        r['subject']['age_min'] = 67
        r['subject']['age_max'] = 67
        r['subject']['age_unit']['id'] = 'UO:0000036'
        r['subject']['age_unit']['label'] = 'year'
        r['subject']['diagnosis'][0]['study_group_description'] = 'PR'
        r['subject']['diagnosis'][0]['disease_diagnosis']['id'] = 'DOID:4450'
        r['subject']['diagnosis'][0]['disease_diagnosis']['label'] = 'renal cell carcinoma'
    elif r['subject']['subject_id'] == '122':
        r['subject']['sex'] = 'female'
        r['subject']['race'] = 'White'
        r['subject']['age_min'] = 41
        r['subject']['age_max'] = 41
        r['subject']['age_unit']['id'] = 'UO:0000036'
        r['subject']['age_unit']['label'] = 'year'
        r['subject']['diagnosis'][0]['study_group_description'] = 'CR'
        r['subject']['diagnosis'][0]['disease_diagnosis']['id'] = 'DOID:4450'
        r['subject']['diagnosis'][0]['disease_diagnosis']['label'] = 'renal cell carcinoma'
    elif r['subject']['subject_id'] == '128':
        r['subject']['sex'] = 'male'
        r['subject']['race'] = 'White'
        r['subject']['age_min'] = 50
        r['subject']['age_max'] = 50
        r['subject']['age_unit']['id'] = 'UO:0000036'
        r['subject']['age_unit']['label'] = 'year'
        r['subject']['diagnosis'][0]['study_group_description'] = 'CR'
        r['subject']['diagnosis'][0]['disease_diagnosis']['id'] = 'DOID:4450'
        r['subject']['diagnosis'][0]['disease_diagnosis']['label'] = 'renal cell carcinoma'
    else:
        print('unknown subject', r['subject'])

    r['subject']['age_event'] = "enrollment"
    r['subject']['species']['id'] = "NCBITAXON:9606"
    r['subject']['species']['label'] = "Homo sapiens"
    del r['subject']['age']

    r['data_processing'][0]['data_processing_id'] = '0dfadecb-9a30-4182-98bd-f07b8437218c-007'
    r['data_processing'][0]['primary_annotation'] = True
    r['data_processing'][0]['software_versions'] = 'igblast-ls5-1.14u1'
    r['data_processing'][0]['germline_database'] = 'VDJServer IMGT 2019.01.23'

    files = []
    fname = ''
    for entry in r['sample']:
        fname = entry['sequencing_files']['filename']
        fname = fname.replace('.fasta','')
        fname = fname + '.igblast.airr.tsv.gz'
        files.append(fname)
    r['data_processing'][0]['data_processing_files'] = [ ','.join(files) ]

airr.write_repertoire(out_file, reps)
