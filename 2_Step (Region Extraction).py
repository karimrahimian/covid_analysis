import pandas as pd
import numpy as np
import collections
import os

nsp_ref = pd.read_excel("data/info/ref_sequence.xlsx",sheet_name='usage')
nsp_ref = nsp_ref.to_numpy()

ref_country = pd.read_excel("clean/AllCountry.xlsx")
ref_country = ref_country.to_numpy()
ref_country = ref_country[:,1]

col_name = ['ID1', 'ID2', 'Date', 'Length', 'Sequence', 'Country', 'Continent', 'GenType']

def findNanIndex(data):
    T = []
    for i,item in enumerate(data):
        result = isinstance(item, str)
        if result==False:
            if np.isnan(item):
                T.append(i)
    return T

for nsp in nsp_ref:
    print (nsp[0])
    maindata = pd.read_csv("clean/total/{}.csv".format(nsp[0]))
    maindata = maindata.to_numpy()

    continents = maindata[:,6]
    continents = collections.Counter(continents)

    for key in continents.keys():
        data = np.empty((continents[key], 8), dtype=np.object)
        indexes = np.where(maindata[:,6]==key)
        if (not isinstance(key,str)):
            if (np.isnan(key)):
                nanIndex = findNanIndex(maindata[:,6])
                nandata = maindata[nanIndex,:]
                df = pd.DataFrame(nandata,columns =col_name)
                df.to_csv("clean/nan/{}.csv".format(nsp[0]))
        else:
            data = maindata[list(indexes),:]
            data = data.reshape(data.shape[1], data.shape[2])
            df = pd.DataFrame(data,columns=col_name)
            try:
                os.mkdir("clean/region/{}".format(key.capitalize()))
            except:
                pass
            df.to_csv('clean/region/{}/{}.csv'.format(key.capitalize(),nsp[0]),index=False)

    country = maindata[:,5]
    country = collections.Counter(country)

    for key in country.keys():
        if (isinstance(key,str)):
            data = np.empty((country[key], 8), dtype=np.object)
            indexes = np.where(maindata[:,5]==key)
            data = maindata[list(indexes),:]
            data = data.reshape(data.shape[1],data.shape[2])
            col_name = ['ID1', 'ID2', 'Date', 'Length', 'Sequence', 'Country', 'Continent', 'GenType']
            df = pd.DataFrame(data,columns=col_name)
            try:
                os.mkdir("clean/region/{}".format(key))
            except:
                pass
            try:
                df.to_csv('clean/region/{}/{}.csv'.format(key.capitalize(), nsp[0]), index=False)
            except Exception as e:
                print(e)


