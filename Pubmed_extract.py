# author: Ivy Wang
# May, 2019
import pickle
import datetime
from Bio import Entrez
from Bio.Entrez import efetch, esearch, esummary
import xmltodict
import csv
import itertools
from collections import OrderedDict
import nltk
from nltk.tag import StanfordNERTagger
from nltk import ne_chunk, pos_tag, word_tokenize
from nltk.tree import Tree
# import ner
# from bs4 import BeautifulSoup
import json
import re
import string
# import sklearn
# import sklearn_crfsuite
# from sklearn_crfsuite import scorers
# from sklearn_crfsuite import metrics
# import pymongo
# from pymongo import MongoClient
import pandas as pd
import collections
import http.client

Entrez.email = 'zw85@georgetown.edu'


def search_for_id(query):
    """
    queries the disease name and fetch matching PubMed IDs from PubMed. change retmax= parameter to change # of results
    :param query:  string of disease name plus other search items
    :return: list of Pubmed IDs that match
    """
    handle = esearch(db='pubmed', term=query, retmax=10000, rettype='uilist', retmode='json')
    xml_format = handle.read()
    xml_json = json.loads(xml_format)

    id_list = xml_json["esearchresult"]["idlist"]
    return id_list


def print_xml(pmid):
    """
    from the result of the query (pubmed ids) to parse out title, abstract, article type, and journal info
    :param pmid: list of PubMed IDs after querying
    :return: list of titles,abstracts, types, and journal names
    """
    handle = efetch(db='pubmed', id=','.join(map(str, pmid)), retmode='xml', rettype='text')
    # try:
    print('entering print_xml')
    doc = handle.read()
    # except http.client.IncompleteRead:
    #     continue
    doc = xmltodict.parse(doc)
    doc = json.loads(json.dumps(doc))
    print('have read the doc')
    d = doc['PubmedArticleSet']["PubmedArticle"]
    titles = []
    types = []
    abstracts = []
    jour_names = []
    for i in d:  # iterate through each article
        # find journal information
        if i["MedlineCitation"]['Article']['Journal']['Title'] is not None:
            jour_name = i["MedlineCitation"]['Article']['Journal']['Title']
            jour_names.append(jour_name)
        else:
            jour_names.append('no journal found')
        # find title information
        t = i["MedlineCitation"]['Article']['ArticleTitle']
        if isinstance(t, str):
            t = i["MedlineCitation"]['Article']['ArticleTitle']
        elif i["MedlineCitation"]['Article']['ArticleTitle'] is None:
            t = "no title"
        else:
            t = i["MedlineCitation"]['Article']['ArticleTitle']['#text']
        titles.append(t)
        if 'Abstract' in i['MedlineCitation']['Article']:
            abstracts.append(i['MedlineCitation']['Article']['Abstract']['AbstractText'])
        else:
            abstracts.append('no abstract')
        # find type of article
        type = i['MedlineCitation']['Article']['PublicationTypeList']['PublicationType']
        if isinstance(type, dict):
            types.append(type['#text'])
            # print(ty['#text'])
        else:
            # print(ty)
            type_stripped = []
            for d in type:
                type_stripped.append(d['#text'])
            type_stripped = ', '.join(type_stripped)
            types.append(type_stripped)
    return titles, abstracts, types, jour_names



# loading pickled drugbank (filtered list) and journal impact factor information (dict from journal to if),
# see directory for pickled files
pickle_drugs = open("list.pickle", "rb")
drugs = pickle.load(pickle_drugs)
pickle_drugs.close()

pickle_if = open("dict.pickle", "rb")
if_dict = pickle.load(pickle_if)
pickle_if.close()

# possible phase seen
phases = ['stage III', 'stage IV', 'advanced', 'metastasis', 'metastatic', 'matastatic malignant',
          'localized', 'intermediate', 'relapsed', 'refractory', 'BRAF', 'BRAF-mutant',
          ' early-stage', 'initial therapy', ' first-line therapy', 'stage III',
          'inoperable', 'neoadjuvant', 'adjuvant',
          'stage II', 'stage I', 'in-transit', 'induction', 'consolidation',
          'systemic', 'low risk', 'high risk', 'distant nodal',
          'initial', 'second-line', 'subsequent', 'first-line', 'additional therapy', 'premenopausal', 'postmenopausal',
          'HER2-Negative',
          'HER2-Positive', 'HER2-targeted', 'resectable', 'unresectable', 'liver only', 'lung only', 'low-risk',
          'high-risk',
          'Previous chemotherapy', 'KRAS', 'NRAS', 'BRAF V600E', 'intensive', 'Previous oxaliplatin-based',
          'Previous irinotecan-based',
          'Previous FOLFOXIRI', 'Previous ﬂuoropyrimidine', 'dMMR/MSI-H',
          'Platinum-based', 'concurrent systemic', 'concurrent', 'induction', 'inductive', 'end-stage', 'end stage'
                                                                                                        'clear-cell']
# end preprocessing/preparation


# build entries of query diseases; count:
query_list = ['acute lymphoblastic leukemia',
             'acute myeloid leukemia', 'acute myeloid leukemia OR acute promyelocytic leukemia OR APL',
             'systemic light chain amyloidosis',
              'anal carcinoma', 'b-cell lymphomas OR b-cell lymphoma', 'b-cell lymphoma OR follicular lymphoma',
              'b-cell lymphoma OR marginal zone lymphoma', 'b-cell lymphoma OR mantle cell lymphoma',
              'b-cell lymphoma OR diffuse large c-cell lymphoma',
              'b-cell lymphoma OR high-grade b-cell lymphoma',
              'b-cell lymphoma OR burkitt lymphoma',
              'b-cell lymphoma OR aids-related b-cell lymphoma',
              'b-cell lymphoma OR lymphoblastic lymphoma',
              'b-cell lymphoma OR post-transplant lymphoproliferative disorder',
              'b-cell lymphoma OR castleman’s disease',
              'basal cell skin cancer',
              'bladder cancer OR urothelial carcinoma OR urothelial cancer', 'bone cancer OR osteosarcoma',
              'bone cancer OR chondrosarcoma',
              'bone cancer OR chordoma', 'bone cancer OR ewing sarcoma',
               'bone cancer OR giant cell tumor of the bone',
              'breast cancer OR breast carcinoma', 'breast cancer OR ductal carcinoma in situ OR DCIS',
              'breast cancer OR inflammatory breast cancer',
              'cervical cancer',
              'chronic lymphocytic leukemia OR small lymphocytic lymphoma OR CLL/SLL OR CLL',
              'chronic myeloid leukemia','colon cancer OR colorectal cancer', 'cutaneous melanoma',
              'dermatofibrosarcoma protuberans',

              'central nervous system cancer OR glioblastoma OR anaplastic glioma',
              'central nervous system cancer OR adult low-grade glioma OR pilocytic and infiltrative supratentorial astrocytoma OR oligodendroglioma',
              'central nervous system cancer OR adult intracranial and spinal ependymoma',
              'central nervous system cancer OR adult medulloblastoma',
              'central nervous system cancer OR primary CNS lymphoma',
              'central nervous system cancer OR primary spinal cord tumor',
              'central nervous system cancer OR meningioma',
              'central nervous system cancer OR limited brain metastases',
              'central nervous system cancer OR extensive brain metastases',
              'central nervous system cancer OR leptomeningeal metastases',
              'central nervous system cancer OR metastatic spine tumor',


              'esophageal cancer OR esophagogastric junction cancer',
              'esophageal cancer OR gastroesophageal junction cancer',
              'gastric cancer', 'gestational trophoblastic neoplasia', 'hairy cell leukemia',

              '(cancer AND (lip OR mucosa)) OR head and neck cancer',
              'head and neck cancer OR （cancer AND oral cavity)',
              'head and neck cancer OR (oropharynx AND cancer) OR oropharyngeal carcinoma',
              'head and neck cancer OR (hypopharynx AND cancer) OR hypopharyngeal carcinoma',
              'head and neck cancer OR (nasopharynx AND cancer) OR nasopharyngeal carcinoma',
              'head and neck cancer OR (glottic larynx AND cancer)',
              'head and neck cancer OR (cancer AND supraglottic larynx)',
              'head and neck cancer OR (cancer AND ethmoid sinus tumor)',
              'head and neck cancer OR (cancer AND maxillary sinus tumor)',
              'head and neck cancer OR (cancer AND salivary gland tumor) OR salivary duct carcinoma',
              'head and neck cancer OR (cancer AND mucosal melanoma)',


              'hepatobiliary cancer OR hepatocellular carcinoma',
              'hepatobiliary cancer OR gallbladder cancer', 'hepatobiliary cancer OR intrahepatic cholangiocarcinoma',
              'hepatobiliary cancer OR extrahepatic cholangiocarcinoma',

              'hodgkin lymphoma OR classic hodgkin lymphoma',
              'hodgkin lymphoma OR nodular lymphocyte-predominant hodgkin lymphoma',
              'renal cell carcinoma OR kidney cancer', 'malignant pleural mesothelioma', 'merkel cell carcinoma',
              'multiple myeloma',
              'myelodysplastic syndrome','myelodysplastic syndrome OR symptomatic anemia',
              'myeloproliferative neoplasm OR myelofibrosis',

              'neuroendocrine AND adrenal gland tumor',
              'neuroendocrine AND gastrointestinal tract, lung, and thymus tumor OR carcinoid tumor',
              'neuroendocrine tumor AND pancreas',
              'neuroendocrine tumor AND unknown primary',
              'neuroendocrine tumor AND (pheochromocytoma OR paraganglioma)',
              'neuroendocrine tumor AND multiple endocrine neoplasia',

              'non-small cell lung cancer',

              'occult primary OR cancer of unknown primary OR CUP',
              'occult primary OR cancer of unknown primary OR CUP OR adenocarcinoma',
              'occult primary OR cancer of unknown primary OR CUP OR squamous cell carcinoma',

              'ovarian cancer OR epithelial ovarian cancer OR fallopian tube cancer OR primary peritoneal cancer',
              'ovarian cancer OR carcinosarcoma OR malignant mixed müllerian tumor',
              'ovarian cancer OR clear cell carcinoma of the ovary',
              'ovarian cancer OR mucinous carcinoma of the ovary',
              'ovarian cancer OR grade 1 endometrioid carcinoma',
              'ovarian cancer OR low-grade serous carcinoma',
              'ovarian cancer OR ovarian borderline epithelial tumor OR low malignant potential',
              'ovarian cancer OR malignant sex cord-stromal tumor',
              'ovarian cancer OR malignant germ cell tumor',

              'pancreatic adenocarcinoma OR pancreatic cancer OR pancreas cancer', 'penile cancer',
              'primary cutaneous lymphoma OR primary cutaneous marginal zone lymphoma',
              'primary cutaneous lymphoma OR primary cutaneous follicle center lymphoma',
              'primary cutaneous lymphoma OR mycosis fungoides OR sezary syndrome',
              'primary cutaneous lymphoma OR primary cutaneous CD30+ T-cell lymphoproliferative disorder',
              'primary cutaneous lymphoma OR primary cutaneous anaplastic large-cell lymphoma',
              'primary cutaneous lymphoma OR lymphomatoid papulosis',

              'prostate cancer', 'rectal cancer',
              'sarcoma OR extremity superficial trunk soft tissue sarcoma OR head neck soft tissue sarcoma',
              'sarcoma OR retroperitoneal intra-abdominal soft tissue sarcoma',
              'sarcoma OR gastrointestinal stromal tumor',
              'sarcoma OR desmoid tumor OR aggressive fibromatosis',
              'sarcoma OR rhabdomyosarcoma',

              'small cell lung cancer',
              'squamous cell skin cancer',
              'systemic mastocytosis',
              't-cell lymphoma OR anaplastic large cell lymphoma',
              't-cell lymphoma OR peripheral t-cell lymphoma',
              't-cell lymphoma OR breast implant-associated ALCL',
              't-cell lymphoma OR t-cell large granular lymphocytic leukemia',
              't-cell lymphoma OR adult t-cell leukemia OR adult t-cell lymphoma',
              't-cell lymphoma OR t-cell prolymphocytic leukemia',
              't-cell lymphoma OR extranodal NK/T-cell lymphoma nasal type',
              't-cell lymphoma OR hepatosplenic gamma-delta T-cell lymphoma',


              'testicular cancer OR testicular tumor', 'testicular cancer OR seminoma',
              'thymoma OR thymic carcinoma',
              'thyroid carcinoma OR thyroid cancer',
              'thyroid carcinoma OR thyroid cancer OR papillary thyroid carcinoma',
              'thyroid carcinoma OR thyroid cancer OR follicular thyroid carcinoma',
              'thyroid carcinoma OR thyroid cancer OR hürthle cell carcinoma',
              'thyroid carcinoma OR thyroid cancer OR medullary thyroid carcinoma',
              'thyroid carcinoma OR thyroid cancer OR anaplastic thyroid carcinoma',

              'uterine neoplasm OR uterine carcinoma OR endometrial carcinoma OR endometrial cancer',
              'uterine neoplasm OR uterine carcinoma OR uterine sarcoma',

              'uveal melanoma',
              'vulvar cancer OR squamous cell carcinoma',
              'waldenström macroglobulinemia OR waldenstrom macroglobulinemia OR lymphoplasmacytic lymphoma']

# start querying
# for parts below: should enter for loop if doing multiple queries at once, otherwise do as is and change
# the variable disease name
# if entering for loop, properly indent everything afterwards.
# disease_name = 'gastric cancer'

# get today's date
now = str(datetime.datetime.now())
today = now[0:10].replace('-','/')
# print(today)


for disease_name in query_list:
    print(f'for {disease_name}:')
    query = str(disease_name) + f' AND ("2018/01/01"[PDAT] : "{today}"[PDAT])'
    tokens_disease = word_tokenize(disease_name)
    file_name = '_'.join(tokens_disease[:7])
    output_file = str(file_name) + '.csv'

    # find pmids
    pmids_list = search_for_id(query)
    print(f'the pmids are {pmids_list}')
    titles, abstracts, types, journals = print_xml(pmids_list) # get lists of title, abstracts, types, and journal names

    # find impact factor from journal name
    if_list = []
    for j in journals:
        j = j.lower()
        if j in if_dict:
            if_list.append(if_dict[j])
        else:
            if_list.append('n/a')
    print(f'the impact factors:{if_list}')
    # clean up corresponding phase
    phase_list = []
    if 'systemic' in disease_name:
        phases.remove('systemic')
        phase_cleaned = phases
    elif 'systemic' not in phases:
        phases.append('systemic')
        phase_cleaned = phases
    else:
        phase_cleaned = phases

    # find drug name and phases/stages from abstracts
    column_drug_list = []
    stop = ['Goat milk', 'Cow milk', 'Cultivated mushroom', 'dioxygen', 'creatin', 'Mustard seed', 'pear', 'pea', 'fica',
            'oxygen', 'creatinine', 'indole', 'date', 'bia', 'gold', 'veal', 'rosin']
    # entering the abstracts to find drug name, phases
    for i, abs_text in enumerate(abstracts):
        abss = json.dumps(abs_text).lower()  # makes abss a string
        abss.translate(str.maketrans('', '', string.punctuation)).lower()
        abss_tok_list = abss.split()
        # finding drug names
        search = [x for i, x in enumerate(drugs) if x in abss and x not in stop]
        search_long = [x for x in search if len(x) > 5]
        dist_search_long = []
        # for word in search_long:
        #     words = word.split()
        #     if 'treatment' in abss_tok_list:
        #         dist_words = abs(abss_tok_list.index(words[0]) - abss_tok_list.index('treatment'))
        #         # print(dist_words)
        #     if dist_words <= 150:
        #         dist_search_long.append(word)
        if abs_text == 'no abstract':
            search_long = [titles[i]]
        # print(noun_ruled)
        column_drug_list.append(list(set(search_long)))

        # get phases
        tit_match_phase = [x for x in phase_cleaned if x.lower() in titles[i].lower()]
        for p in phase_cleaned:
            if abs_text is not None:
                if p in abs_text:
                    tit_match_phase.append(p)
        phase_list.append(list(set(tit_match_phase)))

    print(f'the phases: {phase_list}')
    print(f'the drugs: {column_drug_list}')

    # getting ready to output to file
    links = []

    for id_number in pmids_list:
        links.append([f"https://www.ncbi.nlm.nih.gov/pmc/articles/{id_number}/"])
    disease_name_column = [disease_name] * (len(pmids_list))
    d = [disease_name_column, pmids_list, titles, links, abstracts, phase_list, column_drug_list, types, if_list]
    header = [['disease name'], ['ID'], ['title'], ['links'], ['abstract'], ['phases'], ['drugs'], ['type of resource'],
              ['impact factor']]
    header_export = zip(*header)
    export_data = zip(*d)

    # write to file
    with open(output_file, 'w', encoding="utf-8", newline='') as c_m:
        wr = csv.writer(c_m)
        wr.writerows(header_export)
        wr.writerows(export_data)
