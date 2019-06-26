import csv
import pandas as pd
import collections
import pickle
# reading in file for impact factor
if_dict = {}
if_file = '/Users/ivywang/PycharmProjects/cancer1/if_index.xlsx'
xl = pd.ExcelFile(if_file)
df1 = xl.parse('Sheet1')
# print(df1)
df = pd.read_excel(if_file)
df = df.where(pd.notnull(df), None)
if_dict  = collections.OrderedDict((k.lower(),v) for (k,v) in df.values)
if_dict['ca-a cancer journal for clinicians'] =  244.585


print(if_dict)
if_pickle_out = open('dict.pickle', 'wb')
pickle.dump(if_dict,if_pickle_out)
if_pickle_out.close()




# drugbank reading in and pickle

fields = ['common_name', 'synonyms']
drug_bank = {key: [] for key in fields}

with open('drugbank_vocabulary.csv', 'r') as drug_bank_file:
    reader = csv.reader(drug_bank_file)
    next(reader)
    for line in reader:
        drug_bank['common_name'].append(line[2].lower())
        if line[5] is not None:
            syns = [x.strip() for x in line[5].lower().split('|')]
            drug_bank['synonyms'].append(syns)
flat_syns = [x for sublist in drug_bank['synonyms'] for x in sublist]
flat_syns = list(filter(None, flat_syns))
drugs = drug_bank['common_name'] + flat_syns
drugs = list(filter(None, drugs))
print(len(drugs))
drugs_pickle_out = open('list.pickle', 'wb')
pickle.dump(drugs,drugs_pickle_out)
drugs_pickle_out.close()