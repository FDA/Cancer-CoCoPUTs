import os

import pandas as pd
from cocoputs_pipeline import pipe

#folder containing rnaseq read count files
folder = ##############################

#desired output folder
#make new folder if doesn't exist yet
out_dir = #############################
if not os.path.isdir(out_dir):
    os.mkdir(out_dir)


#generate all_cancer files from scratch
try:
    all_cancer_dic = {}
    for cancer_type in pipe.DIAG_ORIG_DIC:
        fname_list = pipe.cancer_type_to_fnames(pipe.META_DF, pipe.DIAG_ORIG_DIC[cancer_type])
        try:
            assert len(fname_list)>0

            tpm = pipe.fnames_to_table(folder, fname_list, 'tpm').transpose()
            tpm.index = pd.MultiIndex.from_tuples([(cancer_type, ind) for ind in tpm.index.values])

            all_cancer_dic[cancer_type] = tpm

        except AssertionError:
            pass

    all_cancer_df = pd.concat([all_cancer_dic[c_type] for c_type in all_cancer_dic]).transpose()

    all_cancer_df.to_csv('{o}/all_cancer_tpm.tsv'.format(o=out_dir), sep = '\t')

    for name in pipe.COCOP.keys():
        name_df = pipe.tpm_to_cocop(all_cancer_df, name)
        name_df.to_csv('{o}/all_cancer_{n}_usage.tsv'.format(o=out_dir, n=name), sep = '\t')
except:
    print('Failed all_cancer files')

    
#generate all_normal files from scratch
try:
    all_normal_dic = {}
    for normal_type in pipe.ORIG_DIC:
        fname_list = pipe.normal_type_to_fnames(pipe.META_DF, pipe.ORIG_DIC[normal_type])
        try:
            assert len(fname_list)>0

            tpm = pipe.fnames_to_table(folder, fname_list, 'tpm').transpose()
            tpm.index = pd.MultiIndex.from_tuples([(normal_type, ind) for ind in tpm.index.values])

            all_normal_dic[normal_type] = tpm

        except AssertionError:
            pass

    all_normal_df = pd.concat([all_normal_dic[normal_type] for normal_type in all_normal_dic]).transpose()

    all_normal_df.to_csv('{o}/all_normal_tpm.tsv'.format(o=out_dir), sep = '\t')

    for name in pipe.COCOP.keys():
        name_df = pipe.tpm_to_cocop(all_normal_df, name)
        name_df.to_csv('{o}/all_normal_{n}_usage.tsv'.format(o=out_dir, n=name), sep = '\t')
except:
    print('Failed all_normal files')
    
    
#generate paired files from scratch
try:
    cancer_to_normal = {
         'Cholangiocarcinoma':'Bile Duct',
         'Transitional Cell Carcinoma-Bladder':'Bladder',
         'Carcinoma-Breast':'Breast',
         'Duct and Lobular Carcinoma-Breast':'Breast',
         'Ductal Carcinoma-Breast':'Breast',
         'Lobular Carcinoma-Breast':'Breast',
         'Colorectal Adenocarcinoma':'Colon',
         'Left Colorectal Adenocarcinoma':'Left Colon',
         'Right Colorectal Adenocarcinoma':'Right Colon',
         'Adenocarcinoma-Endometrium':'Endometrium',
         'Endometrioid Adenocarcinoma':'Endometrium',
         'Serous cystadenocarcinoma-Endometrium':'Endometrium',
         'Esophageal Adenocarcinoma':'Esophagus',
         'Squamous Cell Carcinoma-Esophagus':'Esophagus',
         'Squamous Cell Carcinoma-Head and Neck':'Head and Neck',
         'Hepatocellular carcinoma':'Liver',
         'Adenocarcinoma with mixed subtypes-Lung':'Lung',
         'Adenocarcinoma-Lung':'Lung',
         'Bronchioloalveolar carcinoma':'Lung',
         'Mucinous Adenocarcinoma-Lung':'Lung',
         'Papillary adenocarcinoma-Lung':'Lung',
         'Squamous Cell Carcinoma-Lung':'Lung',
         'Prostate Adenocarcinoma':'Prostate',
         'Adenocarcinoma-Stomach':'Stomach',
         'Diffuse type carcinoma-Stomach':'Stomach',
         'Intestinal Type Adenocarcinoma-Stomach':'Stomach',
         'Tubular Adenocarcinoma-Stomach':'Stomach',
         'Clear Cell Renal Cell Carcinoma':'Kidney',
         'Chromophobe Renal Cell Carcinoma':'Kidney',
         'Papillary Renal Cell Carcinoma':'Kidney',
         'Renal Cell Carcinoma':'Kidney'
    } 

    paired_dic = {}
    for cancer_type in pipe.DIAG_ORIG_DIC:

        normal_type = cancer_to_normal['-'.join(cancer_type.split(', '))]

        paired_df = pipe.get_paired(pipe.DIAG_ORIG_DIC[cancer_type], pipe.ORIG_DIC[normal_type], pipe.META_DF)
        paired_df['Cancer Type'] = [cancer_type for ind in paired_df.index.values]

        fname_list = list(paired_df.index.values)
        try:
            assert len(fname_list)>0
            paired_tpm = pipe.fnames_to_table(folder, fname_list, 'tpm')
            paired_tpm = paired_tpm.transpose().join(paired_df)
            paired_tpm = paired_tpm.set_index(['Cancer Type', 'Case ID', 'Sample Type'])

            paired_dic[cancer_type] = paired_tpm

        except AssertionError:
            pass

    full_paired_df = pd.concat([paired_dic[c_type] for c_type in paired_dic]).transpose()

    full_paired_df.to_csv('{o}/paired_tpm.tsv'.format(o=out_dir), sep = '\t')

    for name in pipe.COCOP.keys():
        name_df = pipe.tpm_to_cocop(full_paired_df, name)
        name_df.to_csv('{o}/paired_{n}_usage.tsv'.format(o=out_dir, n=name), sep = '\t')
    
except:
    print('Failed paired files')    
    
       
#generate median cancer usage files from scratch
try:
    median_cancer_tpm = {}
    for cancer_type in pipe.DIAG_ORIG_DIC:
        diag_orig = pipe.DIAG_ORIG_DIC[cancer_type]

        cancer_files = pipe.cancer_type_to_fnames(pipe.META_DF, diag_orig)
        median_cancer_tpm[cancer_type] = pipe.fnames_to_median(folder, cancer_files, 'tpm').rename({'median':cancer_type}, axis = 1)

    median_cancer = pd.concat([median_cancer_tpm[x].transpose() for x in median_cancer_tpm]).transpose()

    median_cancer.to_csv('{o}/cancer_median_gene_tpm.tsv'.format(o=out_dir), sep = '\t')

    for name in pipe.COCOP.keys():
        fname = 'cancer_median_gene_{x}_usage.tsv'.format(x=name)
        cancer = pipe.tpm_to_cocop(median_cancer, name).transpose()
        cancer.to_csv('{o}/{f}'.format(o=out_dir, f=fname), sep = '\t')
    
except:
    print('Failed median cancer files')
    
#create median normal usage files from scratch
try:
    median_normal_tpm = {}
    for normal_type in pipe.ORIG_DIC:
        orig = pipe.ORIG_DIC[normal_type]

        normal_files = pipe.normal_type_to_fnames(pipe.META_DF, orig)
        median_normal_tpm[normal_type] = pipe.fnames_to_median(folder, normal_files, 'tpm').rename({'median':normal_type}, axis = 1)

    median_normal = pd.concat([median_normal_tpm[x].transpose() for x in median_normal_tpm]).transpose()

    median_normal.to_csv('{o}/normal_median_gene_tpm.tsv'.format(o=out_dir), sep = '\t')

    for name in pipe.COCOP.keys():
        fname = 'normal_median_gene_{x}_usage.tsv'.format(x=name)
        normal = pipe.tpm_to_cocop(median_normal, name).transpose()
        normal.to_csv('{o}/{f}'.format(o=out_dir, f=fname), sep = '\t')
except:
    print('Failed median normal files')