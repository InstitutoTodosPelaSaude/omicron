## -*- coding: utf-8 -*-

## Created by: Anderson Brito
## Email: anderson.brito@itps.org.br
## Release date: 2022-01-19
## Last update: 2023-03-03
## Refactor by: Bragatte

import pandas as pd
import os
import numpy as np
import hashlib
import time
import argparse
from epiweeks import Week
from tqdm.auto import tqdm



import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)

pd.set_option('display.max_columns', 500)
pd.options.mode.chained_assignment = None

today = time.strftime('%Y-%m-%d', time.gmtime()) ## for snakefile

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Combine and reformat data tables from multiple sources and output a single TSV file",
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

# local run
    # # path = "/Users/**/**/"
    # path = "/Users/*/respat/"
    # input_folder = path + 'data/'
    # rename_file = input_folder + 'rename_columns.xlsx'
    # correction_file = input_folder + 'fix_values.xlsx'
    # cache_file = input_folder + 'combined_cache.tsv'
    # output = input_folder + today + '_combined_data_hlagyn.tsv'


    def load_table(file):
        df = ''
        if str(file).split('.')[-1] == 'tsv':
            separator = '\t'
            df = pd.read_csv(file, encoding='utf-8', sep=separator, dtype='str')
        elif str(file).split('.')[-1] == 'csv':
            separator = ','
            df = pd.read_csv(file, encoding='utf-8', sep=separator, dtype='str')
        elif str(file).split('.')[-1] in ['xls', 'xlsx']:
            df = pd.read_excel(file, index_col=None, header=0, sheet_name=0, dtype='str')
            df.fillna('', inplace=True)
        elif str(file).split('.')[-1] == 'parquet':
            df = pd.read_parquet(file, engine='auto')
            df.fillna('', inplace=True)
        else:
            print('Wrong file format. Compatible file formats: TSV, CSV, XLS, XLSX, PARQUET')
            exit()
        return df
    # print("Done load tables")

    ## load cache file
    if cache_file not in [np.nan, '', None]:
        dfT = load_table(cache_file)
        dfT.fillna('', inplace=True)
    else:
        dfT = pd.DataFrame()
    ## print("Done load cache")

    ## load renaming patterns
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
    # print("Done create dictionary corrections")


    ## load value corrections
    dfC = load_table(correction_file)
    dfC.fillna('', inplace=True)
    dfC = dfC[dfC['lab_id'].isin(["HLAGyn", "any"])] # filter to correct data into fix_values HLAGyn
    # print("Done load corrections")

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
    # print("Load rename_columns")


    def generate_id(column_id):
        id = hashlib.sha1(str(column_id).encode('utf-8')).hexdigest()
        return id
    # print("Done generate hashlib")


    def deduplicate(dfL, dfN, id_columns, test_name):
        ## generate sample id
        dfL['unique_id'] = dfL[id_columns].astype(str).sum(axis=1)  ## combine values in rows as a long string
        dfL['sample_id'] = dfL['unique_id'].apply(lambda x: generate_id(x)[:16])  ## generate alphanumeric sample id with 16 characters

        ## prevent reprocessing of previously processed samples
        if cache_file not in [np.nan, '', None]:
            duplicates = set(dfL[dfL['sample_id'].isin(dfT['sample_id'].tolist())]['sample_id'].tolist())
            if len(duplicates) == len(set(dfL['sample_id'].tolist())):
                print('\n\t\t * ALL samples (%s) were already previously processed. All set!' % test_name)
                dfN = pd.DataFrame()  ## create empty dataframe, and populate it with reformatted data from original lab dataframe
                dfL = pd.DataFrame()

                return dfN, dfL
            else:
                print('\n\t\t * A total of %s out of %s samples (%s) were already previously processed.' % (str(len(duplicates)), str(len(set(dfL['sample_id'].tolist()))), test_name))
                new_samples = len(set(dfL['sample_id'].tolist())) - len(duplicates)
                print('\t\t\t - Processing %s new samples...' % (str(new_samples)))
                dfL = dfL[~dfL['sample_id'].isin(dfT['sample_id'].tolist())]  ## remove duplicates
        else:
            new_samples = len(dfL['sample_id'].tolist())
            print('\n\t\t\t - Processing %s new samples (%s)...' % (str(new_samples), test_name))

        # print('2')
        # print(dfL.head())
        return dfL, dfN
    # print("Done deduplication")
    
    ## Fix datatables
    print('\nFixing datatables...')
    def fix_datatable(dfL,file):
        dfN = dfL
        # print(dfL.columns.tolist())
        # print(''.join(dfL.columns.tolist()))
        # if lab == 'HLAGyn':
        if 'H1N1' in ''.join(dfL.columns.tolist()) or 'Influenza' in ''.join(dfL.columns.tolist()): ## column with unique row data
            test_name = "Painel viral HLAGyn"

            # print('\t\tDados resp_vir >> Correct format. Proceeding...' % test_name)
            # add sample_id and test_kit
            dfL.insert(1, 'sample_id', '')
            dfL.fillna('', inplace=True)

            id_columns = [
                'Pedido',
                'Idade',
                'Sexo',
                'Data Coleta',
                'Cidade',
                'UF',
                ]

            test_columns = []
            if 'Parainfluenza' in ''.join(dfL.columns.tolist()):
                test_columns = [
                    'VIRUS_Influenza A',
                    'VIRUS_Influenza H1N1',
                    'VIRUS_Influenza H3',
                    'VIRUS_Influenza B',
                    'VIRUS_Metapneumovírus',
                    'VIRUS_Sincicial A',
                    'VIRUS_Sincicial B',
                    'VIRUS_Rinovírus',
                    'VIRUS_Parainfluenza 1',
                    'VIRUS_Parainfluenza 2',
                    'VIRUS_Parainfluenza 3',
                    'VIRUS_Parainfluenza 4',
                    'VIRUS_Adenovirus',
                    'VIRUS_Bocavirus',
                    'VIRUS_CoV-229E',
                    'VIRUS_CoV-HKU',
                    'VIRUS_CoV-NL63',
                    'VIRUS_CoV-OC43',
                    'VIRUS_SARS_Like',
                    'VIRUS_SARS-CoV-2',
                    'VIRUS_Enterovírus',
                    'BACTE_Bordetella pertussis',
                    'BACTE_Bordetella parapertussis',
                    'BACTE_Mycoplasma pneumoniae',
                    ]
                dfL.insert(1, 'test_kit', 'test_24')
            elif 'PH4' in ''.join(dfL.columns.tolist()):
                test_columns = [
                    'VIRUS_IA',
                    'VIRUS_H1N1',
                    'VIRUS_AH3',
                    'VIRUS_B',
                    'VIRUS_MH',
                    'VIRUS_SA',
                    'VIRUS_SB',
                    'VIRUS_RH',
                    'VIRUS_PH',
                    'VIRUS_PH2',
                    'VIRUS_PH3',
                    'VIRUS_PH4',
                    'VIRUS_ADE',
                    'VIRUS_BOC',
                    'VIRUS_229E',
                    'VIRUS_HKU',
                    'VIRUS_NL63',
                    'VIRUS_OC43',
                    'VIRUS_SARS',
                    'VIRUS_COV2',
                    'VIRUS_EV',
                    'BACTE_BP',
                    'BACTE_BPAR',
                    'BACTE_MP',
                    ]
                dfL.insert(1, 'test_kit', 'test_24')
            elif 'SARS-CoV-2' in ''.join(dfL.columns.tolist()):
                test_columns = [
                    'Vírus Influenza A',
                    'Vírus Influenza B',
                    'Vírus Sincicial Respiratório A/B',
                    'Coronavírus SARS-CoV-2',
                    ]
                dfL.insert(1, 'test_kit', 'test_4')
            else:
                print('WARNING! Unknown file format. Check for inconsistencies.')

            for column in id_columns + test_columns:
                if column not in dfL.columns.tolist():
                    dfL[column] = ''
                    print('\t\t\t - No \'%s\' column found. Please check for inconsistencies. Meanwhile, an empty \'%s\' column was added.' % (column, column))


            ## assign id and deduplicate
            dfL, dfN = deduplicate(dfL, dfN, id_columns, test_name)
            # print('3')
            # print(dfL.head())

            if dfL.empty:
                # print('# Returning an empty dataframe')
                return dfN


            # starting reformatting process
            pathogens = {
                'SC2': [
                    'VIRUS_SARS_Like',
                    'VIRUS_SARS-CoV-2',
                    'Coronavírus SARS-CoV-2',
                    'VIRUS_SARS',
                    'VIRUS_COV2',
                    ],
                'FLUA': [
                    'VIRUS_Influenza A',
                    'VIRUS_Influenza H1N1',
                    'VIRUS_Influenza H3',
                    'Vírus Influenza A',
                    'VIRUS_IA',
                    'VIRUS_H1N1',
                    'VIRUS_AH3',
                    ],
                'FLUB': [
                    'VIRUS_Influenza B',
                    'Vírus Influenza B',
                    'VIRUS_B',
                    ],
                'VSR': [
                    'VIRUS_Sincicial A',
                    'VIRUS_Sincicial B',
                    'Vírus Sincicial Respiratório A/B',
                    'VIRUS_SA','VIRUS_SB',
                    ],
                'META': [
                    'VIRUS_Metapneumovírus',
                    'VIRUS_MH',
                    ],
                'RINO': [
                    'VIRUS_Rinovírus',
                    'VIRUS_RH',
                    ],
                'PARA': [
                    'VIRUS_Parainfluenza 1',
                    'VIRUS_Parainfluenza 2',
                    'VIRUS_Parainfluenza 3',
                    'VIRUS_Parainfluenza 4',
                    'VIRUS_PH',
                    'VIRUS_PH2',
                    'VIRUS_PH3',
                    'VIRUS_PH4',
                    ],
                'ADENO': [
                    'VIRUS_Adenovirus',
                    'VIRUS_ADE',
                    ], 
                'BOCA': [
                    'VIRUS_Bocavirus',
                    'VIRUS_BOC',
                    ],
                'COVS': [
                    'VIRUS_CoV-229E',
                    'VIRUS_CoV-HKU',
                    'VIRUS_CoV-NL63',
                    'VIRUS_CoV-OC43',
                    'VIRUS_229E',
                    'VIRUS_HKU',
                    'VIRUS_NL63',
                    'VIRUS_OC43',
                    ],
                'ENTERO': [
                    'VIRUS_Enterovírus',
                    'VIRUS_EV',
                    ],
                'BAC': [
                    'BACTE_Bordetella pertussis',
                    'BACTE_Bordetella parapertussis',
                    'BACTE_Mycoplasma pneumoniae',
                    'BACTE_BP',
                    'BACTE_BPAR',
                    'BACTE_MP',
                    ]
                }

            # adding missing columns
            if 'Dt. Nascimento' not in dfL.columns.tolist():
                dfL['birthdate'] = ''
            dfL['Ct_FluA'] = ''
            dfL['Ct_FluB'] = ''
            dfL['Ct_VSR'] = ''
            dfL['Ct_RDRP'] = ''
            dfL['Ct_geneS'] = ''
            dfL['Ct_geneN'] = ''
            dfL['Ct_ORF1ab'] = ''
            dfL['Ct_geneE'] = ''

            dfN = pd.DataFrame()
            unique_cols = list(set(dfL.columns.tolist()))
            for idx, row in dfL.iterrows():
                data = {} # one data row for each request
                for col in unique_cols:
                    value = dfL.loc[idx, col]
                    # data[col] = dfL[col].tolist()[0]
                    data[col] = value

                for pat, tests in pathogens.items():
                    # print('\n' + pat, tests)
                    data[pat + '_test_result'] = 'NA' #'Not tested'

                    for target in tests:
                        # print(file, target)
                        if target in data.keys():
                            # pedido = data['Pedido']
                            # print(pedido, target, data[target])
                            if data[target] == 'Detectado':
                                if data[pat + '_test_result'] in ['NA', 'Não Detectado']:
                                    # print('\t\t\t >>>' + target, data[target], '*', data[pat + '_test_result'])
                                    data[pat + '_test_result'] = 'Pos' #'Detectado'
                                    # print(pedido + ': ' + pat + '_test_result', 'fixed >>>', data[pat + '_test_result'])
                                # else:
                                    # print(pat + '_test_result', 'fixed >>>', data[pat + '_test_result'])
                            elif data[target] == 'Não Detectado':
                                if data[pat + '_test_result'] == 'NA': #'Not tested':
                                    # print('\t\t\t >>>' + target, data[target], '*', data[pat + '_test_result'])
                                    data[pat + '_test_result'] = 'Neg' #'Não Detectado'
                                    # print(pat + '_test_result', 'fixed >>>', data[pat + '_test_result'])
                            else:
                                pass # if neither "Detectado" nor "Não Detectado", live it as is ("NA")

                dfN = dfN.append(data, ignore_index=True)


        elif 'CT_N' in dfL.columns.tolist():
            # print('\t\tDados covid >> Correct format. Proceeding...')
            test_name = "Covid-19"

            id_columns = [
                'Pedido',
                'Idade',
                'Sexo',
                'Data Coleta',
                'Cidade',
                'UF',
                'CT_I',
                'CT_N',
                'CT_ORF1AB',
                ]

            for column in id_columns:
                if column not in dfL.columns.tolist():
                    dfL[column] = ''
                    print('\t\t\t - No \'%s\' column found. Please check for inconsistencies. Meanwhile, an empty \'%s\' column was added.' % (column, column))

            ## generate sample id
            dfL.insert(1, 'sample_id', '')
            dfL.insert(1, 'test_kit', 'covid')
            dfL.fillna('', inplace=True)


            ## assign id and deduplicate
            dfL, dfN = deduplicate(dfL, dfN, id_columns, test_name)
            #print(dfL.head())


            if dfL.empty:
                #print('# Returning an empty dataframe')
                return dfN

            # starting lab specific reformatting
            pathogens = {
                'SC2': [], 
                'FLUA': [],
                'FLUB': [], 
                'VSR': [], 
                'META': [], 
                'RINO': [],
                'PARA': [], 
                'ADENO': [], 
                'BOCA': [], 
                'COVS': [], 
                'ENTERO': [], 
                'BAC': [],
                }

            for p, t in pathogens.items():
                if p != 'SC2':
                    dfL[p + '_test_result'] = 'NA' #'Not tested'

            ## adding missing columns
            if 'Dt. Nascimento' not in dfL.columns.tolist():
                dfL['birthdate'] = ''
            dfL['Ct_FluA'] = ''
            dfL['Ct_FluB'] = ''
            dfL['Ct_VSR'] = ''
            dfL['Ct_RDRP'] = ''
            dfL['Ct_geneS'] = ''
            dfL['Ct_geneN'] = ''
            dfL['Ct_geneE'] = ''

            ## fix decimal values
            dfL['CT_N'] = dfL['CT_N'].str.replace(',', '.')
            dfL['CT_ORF1AB'] = dfL['CT_ORF1AB'].str.replace(',', '.')
            dfN = dfL

        else:
            # print('\t\tFile = ' + file)
            print('\t\tWARNING! Unknown file format. Check for inconsistencies.')
            # exit()
        return dfN
    # print("Done reformating")

    def rename_columns(id, df):
        # print(df.columns.tolist())
        # print(dict_rename[id])
        if id in dict_rename:
            df = df.rename(columns=dict_rename[id])
            # print(df.columns.tolist())
        return df
    # print("Done rename")

    ## fix data points create dict from dict
    def fix_data_points(id, col_name, value):
        new_value = value
        # print(value)
        if value in dict_corrections[id][col_name]:
            new_value = dict_corrections[id][col_name][value]
            # print(str("test = ") + new_value)
        return new_value
    # print("Load dictionary for fix_values")

    ## open data files
    for element in os.listdir(input_folder):
        if not element.startswith('_'):
            if element == 'HLAGyn': ## check if folder is the correct one
                id = element
                element = element + '/'
                if os.path.isdir(input_folder + element) == True:
                    print('\n# Processing datatables from: ' + id)
                    for filename in sorted(os.listdir(input_folder + element)):
                        if filename.split('.')[-1] in ['tsv', 'csv', 'xls', 'xlsx', 'parquet'] and filename[0] not in ['~', '_']:
                            print('\n\t- File: ' + filename)
                            df = load_table(input_folder + element + filename)
                            df.fillna('', inplace=True)
                            df.reset_index(drop=True)

                            df = fix_datatable(df, filename) ## reformat datatable
                            if df.empty:
                                # print('##### Nothing to be done')
                                continue

                            df.insert(0, 'lab_id', id)
                            df = rename_columns(id, df) ## fix data points
                            dfT = dfT.reset_index(drop=True)
                            df = df.reset_index(drop=True)

                            # print(df.head(2)) # only reformat
                            # print(dfT.head(2)) # all labs empty without cache

                            print('\n# Fixing data points...')
                            for lab_id, columns in tqdm(dict_corrections.items()):
                                print('\t- Fixing data from: ' + lab_id)
                                for column, values in columns.items():
                                    # print('\t- ' + column + ' (' + column + ' → ' + str(values) + ')')
                                    df[column] = df[column].apply(lambda x: fix_data_points(lab_id, column, x))
                            
                            # checking duplicates
                            # print(df.columns[df.columns.duplicated(keep=False)])
                            # print(dfT.columns[dfT.columns.duplicated(keep=False)])

                            ## add age from birthdate, if age is missing
                            if 'birthdate' in df.columns.tolist():
                                for idx, row in tqdm(df.iterrows()):
                                    birth = df.loc[idx, 'birthdate']
                                    test = df.loc[idx, 'date_testing']
                                    if birth not in [np.nan, '', None]:
                                        birth = pd.to_datetime(birth)
                                        test = pd.to_datetime(test) ## add to correct dtypes for calculations
                                        age = (test - birth) / np.timedelta64(1, 'Y')
                                        df.loc[idx, 'age'] = np.round(age, 1) ## this gives decimals
                                        #df.loc[idx, 'age'] = int(age)
                                    print(f'Processing tests {idx + 1} of {len(df)}') ## print processed lines 

                                    ## Change the data type of the 'age' column to integer
                                    df['age'] = pd.to_numeric(df['age'], downcast='integer',errors='coerce').fillna(-1).astype(int)
                                    df['age'] = df['age'].apply(int)

                            ## fix sex information
                            df['sex'] = df['sex'].apply(lambda x: x[0] if x != '' else x)

                            frames = [dfT, df]
                            df2 = pd.concat(frames).reset_index(drop=True)
                            dfT = df2

    dfT = dfT.reset_index(drop=True)
    dfT.fillna('', inplace=True)
    # print('Done fix tables')

    ## old place where add `Fixing data points...`

    # reformat dates and get ages
    dfT['date_testing'] = pd.to_datetime(dfT['date_testing'])

    ## create epiweek column
    def get_epiweeks(date):
        try:
            date = pd.to_datetime(date)
            epiweek = str(Week.fromdate(date, system="cdc")) # get epiweeks
            year, week = epiweek[:4], epiweek[-2:]
            epiweek = str(Week(int(year), int(week)).enddate())
#             epiweek = str(Week.fromdate(date, system="cdc"))  # get epiweeks
#             epiweek = epiweek[:4] + '_' + 'EW' + epiweek[-2:]
        except:
            epiweek = ''
        return epiweek

    dfT['epiweek'] = dfT['date_testing'].apply(lambda x: get_epiweeks(x))

    ## Add gene detection results
    def check_detection(ctValue):
        try:
            if ctValue[0].isdigit() and float(ctValue) > 0:
                result = 'Pos' #'Detected'
            elif ctValue[0].isdigit() and float(ctValue) < 1:
                result = 'Neg' #'Not detected'
            else:
                result = 'NA'
        except:
            result = ''
            pass
        return result


    # Ct value columns
    targets = []
    for col in dfT.columns.tolist():
        if col == 'Ct_geneS':
            new_col = col.split('_')[1] + '_detection'
            if new_col not in targets:
                targets.append(new_col)
            dfT[new_col] = dfT[col].apply(lambda x: check_detection(x))

    # reset index
    dfT = dfT.reset_index(drop=True)
    key_cols = [
        'lab_id',
        'test_id',
        'test_kit',
        'sample_id',
        #'region',
        'state',
        'location',
        'date_testing',
        'epiweek',
        'age',
        'sex',
        'FLUA_test_result',
        'Ct_FluA',
        'FLUB_test_result',
        'Ct_FluB',
        'VSR_test_result',
        'Ct_VSR',
        'SC2_test_result',
        'Ct_geneE',
        'Ct_geneN',
        'Ct_geneS',
        'Ct_ORF1ab',
        'Ct_RDRP',
        'geneS_detection',
        'META_test_result',
        'RINO_test_result',
        'PARA_test_result',
        'ADENO_test_result',
        'BOCA_test_result',
        'COVS_test_result',
        'ENTERO_test_result',
        'BAC_test_result',
        ]
    # key_cols = {
    #     'lab_id': 'category',
    #     'test_id': str,
    #     'test_kit': 'category',
    #     'sample_id': str,
    #     # 'region': str,
    #     'state': 'category',
    #     'location': str,
    #     'date_testing': 'datetime64[ns]',
    #     'epiweek': 'datetime64[ns]',
    #     'age': 'int16',
    #     'sex': 'category',
    #     'FLUA_test_result': 'category',
    #     'Ct_FluA': 'float32',
    #     'FLUB_test_result': 'category',
    #     'Ct_FluB': 'float32',
    #     'VSR_test_result': 'category',
    #     'Ct_VSR': 'float32',
    #     'SC2_test_result': 'category',
    #     'Ct_geneE': 'float32',
    #     'Ct_geneN': 'float32',
    #     'Ct_geneS': 'float32',
    #     'Ct_ORF1ab': 'float32',
    #     'Ct_RDRP': 'float32',
    #     'geneS_detection': 'category',
    #     'META_test_result': 'category',
    #     'RINO_test_result': 'category',
    #     'PARA_test_result': 'category',
    #     'ADENO_test_result': 'category',
    #     'BOCA_test_result': 'category',
    #     'COVS_test_result': 'category',
    #     'ENTERO_test_result': 'category',
    #     'BAC_test_result': 'category'
    #     }

    for col in dfT.columns.tolist():
        if col not in key_cols:
            dfT = dfT.drop(columns=[col])

    dfT = dfT[key_cols]

    dfT['date_testing'] = dfT['date_testing'].apply(lambda x: x.strftime('%Y-%m-%d') if pd.notnull(x) else 'XXXXX')

    ## fix test results with empty data
    # for p in pathogens.keys():
    #     dfT[p + '_test_result'] = dfT[p + '_test_result'].apply(lambda x: 'Negative' if x not in ['Negative', 'Positive'] else x)

    ## output duplicates rows
    duplicates = dfT.duplicated().sum()
    if duplicates > 0:
        mask = dfT.duplicated(keep=False) ## find duplicates
        dfD = dfT[mask]
        output2 = input_folder + 'duplicates.tsv'
        dfD.to_csv(output2, sep='\t', index=False)
        print('\nWARNING!\nFile with %s duplicate entries saved in:\n%s' % (str(duplicates), output2))

    ## drop duplicates
    dfT = dfT.drop_duplicates(keep='last')

    ## sorting by date
    dfT = dfT.sort_values(by=['lab_id', 'test_id', 'date_testing'])


    ## time controller for optimization of functions `def`
    ## example with tqdm for def
    # start = time.time()
    # for load_table in tqdm(range(1), desc='Execution Time'):
    #     load_table
    # end = time.time()
    # print("Execution time for load_table: ", end - start)

    start = time.time()
    load_table
    end = time.time()
    print("Execution time for load_table: ", end - start)

    # start = time.time()
    # generate_id
    # end = time.time()
    # print("Execution time for generate_id: ", end - start)

    # start = time.time()
    # deduplicate
    # end = time.time()
    # print("Execution time for deduplicate: ", end - start)

    # start = time.time()
    # fix_datatable
    # end = time.time()
    # print("Execution time for fix_datatable: ", end - start)

    # start = time.time()
    # rename_columns
    # end = time.time()
    # print("Execution time for rename_columns: ", end - start)

    # start = time.time()
    # fix_data_points
    # end = time.time()
    # print("Execution time for fix_data_points: ", end - start)

    # start = time.time()
    # get_epiweeks
    # end = time.time()
    # print("Execution time for get_epiweeks: ", end - start)

    # output combined dataframe
    dfT.to_csv(output, sep='\t', index=False)
    print('\nData successfully aggregated and saved in:\n%s\n' % output)
