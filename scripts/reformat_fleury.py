# -*- coding: utf-8 -*-

# Created by: Anderson Brito
# Email: anderson.brito@itps.org.br
# Release date: 2022-XX-XX
# Last update: 2022-XX-XX
# Refactor by: Bragatte

import pandas as pd
import os
import numpy as np
import hashlib
import time
import argparse
from epiweeks import Week

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)

pd.set_option('display.max_columns', 500)
pd.options.mode.chained_assignment = None

today = time.strftime('%Y-%m-%d', time.gmtime())

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Reformat and combine data tables from Fleury",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--datadir", required=True, help="Name of the folder containing independent folders for each lab")
    parser.add_argument("--rename", required=False, help="TSV, CSV, or excel file containing new standards for column names")
    parser.add_argument("--correction", required=False, help="TSV, CSV, or excel file containing data points requiring corrections")
    parser.add_argument("--cache", required=False, help="Previously processed data files")
    parser.add_argument("--output", required=True, help="TSV file aggregating all columns listed in the 'rename file'")
    args = parser.parse_args()

    path = os.path.abspath(os.getcwd())
    input_folder = path + '/' + args.datadir + '/'
    rename_file = args.rename
    correction_file = args.correction
    cache_file = args.cache
    output = args.output

    # # path = "/Users/anderson/google_drive/ITpS/projetos_itps/resp_pathogens/analyses/dev/20230428_fleury/"
    # path = "/Users/Anderson/Library/CloudStorage/GoogleDrive-anderson.brito@itps.org.br/Outros computadores/My Mac mini/google_drive/ITpS/projetos_itps/resp_pathogens/analyses/dev/20230428_fleury/"
    # input_folder = path + 'data/'
    # rename_file = input_folder + 'rename_columns.xlsx'
    # correction_file = input_folder + 'fix_values.xlsx'
    # cache_file = '' # path + 'data/cache.tsv'#input_folder + '2022-08-02_combined_data_dasa.tsv'
    # output = input_folder + today + '_combined_data_fleury.tsv'

    def load_table(file):
        df = ''
        if str(file).split('.')[-1] == 'tsv':
            separator = '\t'
            df = pd.read_csv(file, encoding='latin-1', sep=separator, dtype='str')
        elif str(file).split('.')[-1] == 'csv':
            separator = ','
            df = pd.read_csv(file, encoding='latin-1', sep=separator, dtype='str')
        elif str(file).split('.')[-1] in ['xls', 'xlsx']:
            df = pd.read_excel(file, index_col=None, header=0, sheet_name=0, dtype='str')
            df.fillna('', inplace=True)
        else:
            print('Wrong file format. Compatible file formats: TSV, CSV, XLS, XLSX')
            exit()
        return df

    # load cache file
    if cache_file not in [np.nan, '', None]:
        dfT = load_table(cache_file)
        dfT.fillna('', inplace=True)
    else:
        # dfP = pd.DataFrame()
        dfT = pd.DataFrame()

    # load renaming patterns
    dfR = load_table(rename_file)
    dfR.fillna('', inplace=True)


    dict_rename = {}
    # dict_corrections = {}
    for idx, row in dfR.iterrows():
        id = dfR.loc[idx, 'lab_id']
        if id not in dict_rename:
            dict_rename[id] = {}
        old_colname = dfR.loc[idx, 'column_name']
        new_colname = dfR.loc[idx, 'new_name']
        rename_entry = {old_colname: new_colname}
        dict_rename[id].update(rename_entry)

    # load value corrections
    dfC = load_table(correction_file)
    dfC.fillna('', inplace=True)

    dict_corrections = {}
    all_ids = list(set(dfC['lab_id'].tolist()))
    for idx, row in dfC.iterrows():
        lab_id = dfC.loc[idx, 'lab_id']
        colname = dfC.loc[idx, 'column_name']

        old_data = dfC.loc[idx, 'old_data']
        new_data = dfC.loc[idx, 'new_data']
        if old_data + new_data not in ['']:
            labs = []
            if colname == 'any':
                labs = all_ids
            else:
                labs = [lab_id]
            for id in labs:
                if id not in dict_corrections:
                    dict_corrections[id] = {}
                if colname not in dict_corrections[id]:
                    dict_corrections[id][colname] = {}
                data_entry = {old_data: new_data}
                dict_corrections[id][colname].update(data_entry)

    def generate_id(column_id):
        id = hashlib.sha1(str(column_id).encode('utf-8')).hexdigest()
        return id


    def deduplicate(dfL, dfN, id_columns):
        # generate sample id
        dfL['unique_id'] = dfL[id_columns].astype(str).sum(axis=1)  # combine values in rows as a long string
        dfL['sample_id'] = dfL['unique_id'].apply(lambda x: generate_id(x))  # generate alphanumeric sample id

        # prevent reprocessing of previously processed samples
        if cache_file not in [np.nan, '', None]:
            duplicates = set(dfL[dfL['sample_id'].isin(dfT['sample_id'].tolist())]['sample_id'].tolist())
            if len(duplicates) == len(set(dfL['sample_id'].tolist())):
                print('\n\t\t * ALL samples (%s) were already previously processed. All set!' % len(duplicates))
                dfN = pd.DataFrame()  # create empty dataframe, and populate it with reformatted data from original lab dataframe
                dfL = pd.DataFrame()

                # print('1')
                # print(dfL.head())

                return dfN, dfL
            else:
                print('\n\t\t * A total of %s out of %s samples (%s) were already previously processed.' % (str(len(duplicates)), str(len(set(dfL['sample_id'].tolist()))), test_name))
                new_samples = len(set(dfL['sample_id'].tolist())) - len(duplicates)
                print('\t\t\t - Processing %s new samples...' % (str(new_samples)))
                dfL = dfL[~dfL['sample_id'].isin(dfT['sample_id'].tolist())]  # remove duplicates
        else:
            new_samples = len(dfL['sample_id'].tolist())
            print('\n\t\t\t - Processing %s new samples...' % (str(new_samples)))

        # print('2')
        # print(dfL.head())
        return dfL, dfN



    # Fix datatables
    print('\nFixing datatables...')

    pathogens = {'FLUA': ['Covidflursvgx - Influenza A', 'Virusmol, Influenza A', 'Virusmol, Influenza A/H1', 'Virusmol, Influenza A/H1-2009', 'Virusmol, Influenza A/H3', 'Vírus Influenza A (Sazonal)', 'Vï¿½rus Influenza A (Sazonal)', 'Vírus respiratórios - Influenzavirus A'],
                 'FLUB': ['Covidflursvgx - Influenza B', 'Virusmol, Influenza B', 'Vírus respiratórios - Influenzavirus B'],
                 'VSR': ['Covidflursvgx - Vírus Sincicial Respiratório', 'Covidflursvgx - Vï¿½rus Sincicial Respiratï¿½rio', 'Virusmol, Vírus Sincicial Respiratório', 'Virusmol, Vï¿½rus Sincicial Respiratï¿½rio', 'Vírus Sincial Respiratório', 'Vírus respiratório - Sincicial', 'Vï¿½rus respiratï¿½rio - Sincicial', 'Vírus respiratórios - Vírus Sincicial Respira'],
                 'SC2': ['Covid 19, Antígeno, teste rápido', 'Covid 19, Antï¿½geno, teste rï¿½pido', 'Covid 19, Detecção por PCR', 'Covid 19, Detecï¿½ï¿½o por PCR', 'Covidflursvgx - SARS-CoV-2', 'Virusmol, SARS-CoV-2'],
                 'META': ['Virusmol, Metapneumovírus Humano', 'Virusmol, Metapneumovï¿½rus Humano'],
                 'PARA': ['Virusmol, Parainfluenza 1', 'Virusmol, Parainfluenza 2', 'Virusmol, Parainfluenza 3', 'Virusmol, Parainfluenza 4', 'Vírus respiratórios - Parainfluenzavirus 1', 'Vírus respiratórios - Parainfluenzavirus 2', 'Vírus respiratórios - Parainfluenzavirus 3', 'Vírus respiratórios - Parainfluenzavirus 4'],
                 'ADENO': ['Virusmol, Adenovírus', 'Virusmol, Adenovï¿½rus', 'Vírus respiratórios - Adenovírus'],
                 'COVS': ['Virusmol, Coronavírus 229E', 'Virusmol, Coronavï¿½rus 229E', 'Virusmol, Coronavírus HKU1', 'Virusmol, Coronavï¿½rus HKU1', 'Virusmol, Coronavírus NL63', 'Virusmol, Coronavï¿½rus NL63', 'Virusmol, Coronavírus OC43', 'Virusmol, Coronavï¿½rus OC43'],
                 'RINO': ['Virusmol, Rinovï¿½rus/Enterovï¿½rus', 'Virusmol, Rinovírus/Enterovírus'],
                 'ENTERO': [],
                 'BOCA': [],
                 'BAC': ['Virusmol, Bordetella parapertussis', 'Virusmol, Bordetella pertussis', 'Virusmol, Chlamydophila pneumoniae', 'Virusmol, Mycoplasma pneumoniae']}

    def fix_datatable(dfL, pathogens):
        # ignore = ['Influenza A e B - teste rápido', 'Virusmol, Rinovírus/Enterovírus', '']

        # remove useless lines from panel tests AGRESPVIR
        ignore = ['Vï¿½rus respiratï¿½rios - detecï¿½ï¿½o', 'INCONCLUSIVO', '']
        dfL = dfL[~dfL['PATOGENO'].isin(ignore)]
        dfL = dfL[~dfL['RESULTADO'].isin(ignore)]


        # generate sample id
        dfL.insert(1, 'sample_id', '')
        dfL.fillna('', inplace=True)

        # dfN = dfL
        dfN = pd.DataFrame()  # create empty dataframe, and populate it with reformatted data from original lab dataframe
        id_columns = ['CODIGO REQUISICAO', 'PACIENTE', 'IDADE', 'SEXO', 'DATA COLETA', 'MUNICIPIO', 'ESTADO']

        for column in id_columns:
            if column not in dfL.columns.tolist():
                dfL[column] = ''
                print(
                    '\t\t\t - No \'%s\' column found. Please check for inconsistencies. Meanwhile, an empty \'%s\' column was added.' % (column, column))

        # missing columns
        # adding missing columns
        dfL['birthdate'] = ''
        dfL['Ct_FluA'] = ''
        dfL['Ct_FluB'] = ''
        dfL['Ct_VSR'] = ''
        dfL['Ct_RDRP'] = ''
        dfL['Ct_geneE'] = ''
        dfL['Ct_ORF1ab'] = ''
        dfL['Ct_geneN'] = ''
        dfL['Ct_geneS'] = ''
        dfL['geneS_detection'] = ''

        # assign id and deduplicate
        dfL, dfN = deduplicate(dfL, dfN, id_columns)

        if dfL.empty:
            return dfN

        # starting lab specific reformatting

        alltypes = [] # DELETE!!!!
        positives = ['DETECTADO (POSITIVO)', 'Influenza A - POSITIVO', 'Influenza B - POSITIVO', 'P O S I T I V O', 'POSITIVO'] # many ways to report positive results

        for i, (code, dfR) in enumerate(dfL.groupby('CODIGO REQUISICAO')):
            data = {} # one data row for each request
            target_pathogen = {}
            for p, t in pathogens.items(): # set all tests as 'NA' first, to than changed it to Pos or Neg as appropriate
                data[p + '_test_result'] = 'NA'
                for g in t:
                    target_pathogen[g] = p

            exam = list(set(dfR['EXAME'].tolist()))[0] # exam name
            if exam not in alltypes:
                alltypes.append(exam)

            # print('\n\n', code, exam, dfR['PATOGENO'].tolist())

            # iterate over different test types
            if exam in ['2019NCOV', 'AGCOVIDNS', 'COVID19GX', 'COVID19POCT', 'COVID19SALI','INFLUENZAPCR', 'AGSINCURG', 'VRSAG']: # single tests
                single_tests = {'2019NCOV': 'covid_pcr', 'AGCOVIDNS': 'covid_antigen', 'COVID19POCT': 'covid_pcr',
                                'COVID19SALI': 'covid_pcr', 'COVID19GX': 'covid_pcr', 'INFLUENZAPCR': 'flu_pcr',
                                'AGSINCURG': 'vsr_antigen', 'VRSAG': 'vsr_antigen'}

                dfR.insert(1, 'test_kit', single_tests[exam])
                dfR['pathogen'] = dfR['PATOGENO'].apply(lambda x: target_pathogen[x])
                pathogen = list(set(dfR['pathogen'].tolist()))[0]

                # assign test result
                test_result = list(set(dfR['RESULTADO'].tolist()))[0]
                if test_result in positives:
                    result = 'Pos'
                    data[pathogen + '_test_result'] = result
                else: # target not detected
                    result = 'Neg'
                    data[pathogen + '_test_result'] = result


            elif exam in ['AGINFLU']: # Influenza viral panel
                dfR.insert(1, 'test_kit', 'flu_antigen')
                dfR = pd.concat([dfR]*2)
                dfR = dfR.reset_index(drop=True)
                dfR.loc[0, 'PATOGENO'] = 'FLUA'
                dfR.loc[1, 'PATOGENO'] = 'FLUB'
                dfR = dfR.rename(columns={'PATOGENO': 'pathogen'})

                if 'Influenza A e B - POSITIVO' in list(set(dfR['RESULTADO'].tolist())): # both are positive
                    dfR.loc[0, 'RESULTADO'] = 'Influenza A - POSITIVO'
                    dfR.loc[1, 'RESULTADO'] = 'Influenza B - POSITIVO'
                elif 'Influenza A - POSITIVO' in list(set(dfR['RESULTADO'].tolist())): # A positive, B negative
                    dfR.loc[0, 'RESULTADO'] = 'Influenza A - POSITIVO'
                    dfR.loc[1, 'RESULTADO'] = 'NEGATIVO'
                elif 'Influenza B - POSITIVO' in list(set(dfR['RESULTADO'].tolist())): # B positive, A negative
                    dfR.loc[0, 'RESULTADO'] = 'NEGATIVO'
                    dfR.loc[1, 'RESULTADO'] = 'Influenza B - POSITIVO'

                # assign test result
                for pathogen, dfG in dfR.groupby('pathogen'):
                    test_result = list(set(dfG['RESULTADO'].tolist()))[0]

                    if test_result in positives:
                        result = 'Pos'
                        data[pathogen + '_test_result'] = result
                    else:  # target not detected
                        result = 'Neg'
                        data[pathogen + '_test_result'] = result


            elif exam in ['COVIDFLURSVGX']: # viral panel with 4 targets
                dfR.insert(1, 'test_kit', 'test_4')
                dfR['pathogen'] = dfR['PATOGENO'].apply(lambda x: target_pathogen[x])

                for pathogen, dfG in dfR.groupby('pathogen'):
                    test_result = list(set(dfG['RESULTADO'].tolist()))[0]

                    if test_result in positives:
                        result = 'Pos'
                        data[pathogen + '_test_result'] = result
                    else:  # target not detected
                        result = 'Neg'
                        data[pathogen + '_test_result'] = result


            elif exam in ['AGRESPVIR']:  # viral panel with 7 targets
                pass # this block still needs to be developed


            elif exam in ['VIRUSMOL']:  # panel 20 pathogens
                dfR.insert(1, 'test_kit', 'test_21')
                dfR['pathogen'] = dfR['PATOGENO'].apply(lambda x: target_pathogen[x])

                # assign test result
                for pathogen, dfG in dfR.groupby('pathogen'):
                    # print(dfG[['pathogen', 'RESULTADO']])
                    test_results = dfG['RESULTADO'].tolist()

                    positive = ''
                    for res in test_results:
                        if positive == '':
                            # print(pathogen, res)
                            if res in positives:
                                result = 'Pos'
                                data[pathogen + '_test_result'] = result
                                positive = 'yes'
                            else:  # target not detected
                                result = 'Neg'
                                data[pathogen + '_test_result'] = result
            else:
                pass
                # print('Unknown test = ' + exame)

            # add remaining columns into data dict
            unique_cols = list(set(dfR.columns.tolist()))
            # print(unique_cols)

            for col in unique_cols:
                data[col] = dfR[col].tolist()[0]
            # print(data)

            dfN = dfN.append(data, ignore_index=True)
            # print(set(dfN['COVS_test_result'].tolist()))

        # dfN = dfN[['CODIGO REQUISICAO', 'EXAME', 'PATOGENO', 'RESULTADO', 'sample_id', 'FLUA_test_result', 'FLUB_test_result', 'VSR_test_result', 'SC2_test_result', 'META_test_result', 'RINO_test_result', 'PARA_test_result', 'ADENO_test_result', 'BOCA_test_result', 'COVS_test_result', 'ENTERO_test_result', 'BAC_test_result']]
        # dfN = dfN.sort_values(by=['EXAME'])
        # dfN.to_csv(output, sep='\t', index=False)
        # print(alltypes)
        return dfN


    def rename_columns(id, df):
        # print(df.columns.tolist())
        # print(dict_rename[id])
        if id in dict_rename:
            df = df.rename(columns=dict_rename[id])
        # print(df.columns.tolist())
        return df

    # open data files
    for element in os.listdir(input_folder):
        if not element.startswith('_'):
            if element == 'Fleury': # check if folder is the correct one
                id = element
                element = element + '/'
                if os.path.isdir(input_folder + element) == True:
                    print('\n# Processing datatables from: ' + id)
                    for filename in sorted(os.listdir(input_folder + element)):
                        if filename.split('.')[-1] in ['tsv', 'csv', 'xls', 'xlsx'] and filename[0] not in ['~', '_']:
                            print('\n\t- File: ' + filename)
                            df = load_table(input_folder + element + filename)
                            df.fillna('', inplace=True)
                            df.reset_index(drop=True)

                            # print('- ' + str(len(dfT['sample_id'].tolist())) + ' samples (pre)')
                            # print('- ' + str(len(dfT['sample_id'].tolist())) + ' samples (pre)')

                            df = fix_datatable(df, pathogens) # reformat datatable

                            if df.empty:
                                # print('##### Nothing to be done')
                                continue

                            df.insert(0, 'lab_id', id)
                            df = rename_columns(id, df) # fix data points

                            dfT = dfT.reset_index(drop=True)
                            df = df.reset_index(drop=True)

                            frames = [dfT, df]
                            df2 = pd.concat(frames).reset_index(drop=True)
                            dfT = df2
                            # print('##### Concatenating some data')
                            # print('- ' + str(len(dfT['sample_id'].tolist())) + ' samples (pos)')

                            # dfT.to_csv(output, sep='\t', index=False)


    dfT = dfT.reset_index(drop=True)
    dfT.fillna('', inplace=True)

    # fix data points
    def fix_data_points(id, col_name, value):
        new_value = value
        if value in dict_corrections[id][col_name]:
            new_value = dict_corrections[id][col_name][value]
        return new_value

    # print(dfT.head())
    print('\n# Fixing data points...')
    for lab_id, columns in dict_corrections.items():
        print('\t- Fixing data from: ' + lab_id)
        for column, values in columns.items():
            # print('\t- ' + column + ' (' + column + ' → ' + str(values) + ')')
            dfT[column] = dfT[column].apply(lambda x: fix_data_points(lab_id, column, x))

    # reformat dates and convert to datetime format
    dfT['date_testing'] = pd.to_datetime(dfT['date_testing'], dayfirst=True) # , format='%Y-%m-%d', errors='ignore'


    # create epiweek column
    def get_epiweeks(date):
        try:
            date = pd.to_datetime(date)
            epiweek = str(Week.fromdate(date, system="cdc")) # get epiweeks
            year, week = epiweek[:4], epiweek[-2:]
            epiweek = str(Week(int(year), int(week)).enddate())
        except:
            epiweek = ''
        return epiweek

    dfT['epiweek'] = dfT['date_testing'].apply(lambda x: get_epiweeks(x))

    # print(dfT.columns.tolist())
    # print(dfT.head())
    # print(len(dfT['sample_id'].tolist()))

    # add age from birthdate, if age is missing
    if 'birthdate' in dfT.columns.tolist():
        for idx, row in dfT.iterrows():
            birth = dfT.loc[idx, 'birthdate']
            test = dfT.loc[idx, 'date_testing']
            if birth not in [np.nan, '', None]:
                birth = pd.to_datetime(birth)
                age = (test - birth) / np.timedelta64(1, 'Y')
                dfT.loc[idx, 'age'] = np.round(age, 1)

    # fix age
    dfT['age'] = dfT['age'].apply(lambda x: '0' if 'D' in x and 'A' not in x else x)
    dfT['age'] = dfT['age'].apply(lambda x: x.split('A')[0] if 'A' in x else x)

    # fix sex information
    dfT['sex'] = dfT['sex'].apply(lambda x: x[0] if x != '' else x)


    # reset index
    dfT = dfT.reset_index(drop=True)
    # dfT.to_csv(output, sep='\t', index=False)

    key_cols = ['lab_id', 'test_id', 'patient_id', 'test_kit', 'sample_id', 'state', 'location', 'sex', 'date_testing', 'epiweek',
                'age','FLUA_test_result', 'Ct_FluA', 'FLUB_test_result', 'Ct_FluB', 'VSR_test_result', 'Ct_VSR','SC2_test_result', 'Ct_geneE',
                'Ct_geneN', 'Ct_geneS', 'Ct_ORF1ab', 'Ct_RDRP', 'geneS_detection', 'META_test_result', 'RINO_test_result',
                'PARA_test_result', 'ADENO_test_result', 'BOCA_test_result', 'COVS_test_result', 'ENTERO_test_result', 'BAC_test_result']

    for col in dfT.columns.tolist():
        if col not in key_cols:
            dfT = dfT.drop(columns=[col])

    # print(dfT.columns.tolist)

    # keep only key columns, and find null dates
    dfT = dfT[key_cols]
    dfT['date_testing'] = dfT['date_testing'].apply(lambda x: x.strftime('%Y-%m-%d') if pd.notnull(x) else 'XXXXX')

    # # fix test results with empty data
    # for p in pathogens.keys():
    #     dfT[p + '_test_result'] = dfT[p + '_test_result'].apply(lambda x: 'Neg' if x not in ['Neg', 'Pos'] else x)



    # output duplicates rows
    duplicates = dfT.duplicated().sum()
    if duplicates > 0:
        mask = dfT.duplicated(keep=False) # find duplicates
        dfD = dfT[mask]
        output2 = input_folder + 'duplicates.tsv'
        dfD.to_csv(output2, sep='\t', index=False)
        print('\nWARNING!\nFile with %s duplicate entries saved in:\n%s' % (str(duplicates), output2))


    # drop duplicates
    dfT = dfT.drop_duplicates(keep='last')

    # sorting by date
    dfT = dfT.sort_values(by=['lab_id', 'test_id', 'date_testing'])

    # output combined dataframe
    dfT.to_csv(output, sep='\t', index=False)
    print('\nData successfully aggregated and saved in:\n%s\n' % output)
