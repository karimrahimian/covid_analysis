import numpy as np
import pandas as pd
import glob
import os
import pycountry_convert as pc

nsp_ref = pd.read_excel("data/info/ref_sequence.xlsx",sheet_name='usage')
nsp_ref = nsp_ref.to_numpy()

country_ref = pd.read_excel("clean/country_names_final.xlsx")
country_names = country_ref.to_numpy()[:,1]
nsp_names = nsp_ref[:,0]

prefix_path = 'g:/covid_result'
def GetFolders():
    path = glob.glob('clean/region/*')
    dirs = []
    for dir in path:
        dirs.append(dir+"/")
    return dirs

folder_path = GetFolders()

region_statistic = np.zeros((len(country_names),len(nsp_names)))

def Calculate_Statistic():
    for nsp in nsp_ref:
        print(nsp[0])
        for dir in folder_path:
            try:
                maindata = pd.read_csv("{}{}.csv".format(dir,nsp[0]))
                maindata = maindata.to_numpy()
                country_name = dir.split("\\")[1].replace("/","")

                #country_names.add(country_name)

                country_index = np.where(country_names==country_name)
                nsp_index = np.where(nsp_names==nsp[0])
                region_statistic[country_index,nsp_index] = maindata.shape[0]

            except:
                print ("Error on " + dir)
                pass


def country_to_continent(country_name):
    country_alpha2 = pc.country_name_to_country_alpha2(country_name)
    country_continent_code = pc.country_alpha2_to_continent_code(country_alpha2)
    country_continent_name = pc.convert_continent_code_to_continent_name(country_continent_code)
    return country_continent_name

def load_country():
    df_region = pd.read_excel("statistic/statistic_region.xlsx")
    data = df_region.to_numpy()
    col_names = df_region.columns
    col_names = np.insert(col_names, 1, "continents")
    df = pd.read_excel("clean/AllCountry.xlsx")
    cont_ref = df.to_numpy()

    countries = data[:,0]
    continents = np.empty(len(countries),dtype=np.object)

    for i,country in enumerate(countries):
        try:
            index = np.where(cont_ref[:,1]==country)
            continent_name = cont_ref[index[0],2][0]
            continents[i] = cont_ref[index,2][0][0]
        except:
            pass
    data = np.insert(data,0,continents,axis=1)
    df = pd.DataFrame(data,columns=col_names)
    df.to_excel("statistic/statistic_region_country.xlsx")
    print (continents)

#Calculate_Statistic()
#df = pd.DataFrame(region_statistic,index=country_names,columns=nsp_names)
#df.to_excel("statistic/statistic_region.xlsx")

load_country()