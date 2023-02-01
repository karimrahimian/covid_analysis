import pandas as pd
import numpy as np
import collections
import datetime
import json
from dateutil import parser
from sklearn.decomposition import PCA
from collections import OrderedDict
from sklearn.preprocessing import StandardScaler
import seaborn as sns
from matplotlib import pyplot as plt
import collections
import os
import glob


nsp_ref = pd.read_excel("data/info/ref_sequence.xlsx",sheet_name='usage')
nsp_ref = nsp_ref.to_numpy()

aminoacids = np.array(["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K","M", "F", "P", "S", "T", "W", "Y", "V"])
sequence_col_index = 4
folder_path = ["clean/World Wide/",
               "clean/africa/",
               "clean/asia/",
               "clean/europe/",
               "clean/north america/",
               "clean/oceania/",
               "clean/south america/"]

folder_path = ["clean/group_country/NorthAsia/",
               "clean/group_country/WestAsia/"]



#folder_path = ["total/continent/north america/","total/continent/oceania/","total/continent/south america/"]
#sheetname = ["North America","Occeania","South America"]
#sheetname = ["World Wide","Africa","Asia","Europe","North America","Occeania","South America"]
prefix_path = 'g:/covid_result'

def GetFolders():
    path = glob.glob('clean/region/*')
    dirs = []
    for dir in path:
        dirs.append(dir+"/")
    return dirs

folder_path = GetFolders()

def DistictMutant(ref_seq,seq):
    i=0
    mutanttype = ""
    founded =False

    for seq1, seq2 in zip(ref_seq,seq):
        if (seq1 != seq2):
            mutanttype = mutanttype + "{}({})>{},".format(str(seq1).upper(),i+1,str(seq2).upper())
            founded = True
        i+=1
    return founded,mutanttype

def FindMutant(sequence1,sequence2):
    mutants = np.zeros(len(sequence1))
    mutantfound=False
    newmutant = np.zeros((len(sequence2), len(aminoacids)))
    i=0
    for seq1, seq2 in zip(sequence1, sequence2):
        if (seq1 == seq2):
            mutants[i] = 0
        else:
            mutants[i] = 1
            index =np.where( aminoacids == str(seq2).upper())
            newmutant[i,index] = 1
            mutantfound=True
        i+=1
    return mutantfound,mutants ,newmutant

def ConvertToDateTime(date):
    d = date.split('-')
    y = d[0]
    m = d[1]
    d = d[2]
    if (m=="00"):
        m="01"
    if (d=="00"):
        d="01"
    newtime = y+"-"+m+"-"+d
    return datetime.datetime.strptime(newtime, '%Y-%m-%d')

def Find_Heatmap_Data(mutant):
    heat_count = np.zeros(10)
    heat_label = []
    divide_rules=[]
    windowlen = len(mutant) // 10
    residue = len(mutant) % 10

    for i in range(9):
        divide_rules.append(windowlen)

    if (residue>0):
        divide_rules.append(windowlen+residue)
    else:
        divide_rules.append(windowlen)

    #split_mutant = np.split(mutant,divide_rules)
    for index,array_list in enumerate(divide_rules):
        #heat_count[index]= np.sum(array_list)
        if (index==0):
            window = mutant[0:windowlen]
            heat_count[index] = np.sum(window)
            heat_label.append("{}_{}".format(1,windowlen))

        else:
            heat_label.append("{}_{}".format(windowlen * index, (index+1)* windowlen))
            window = mutant[windowlen * index:(index+1)* windowlen]
            heat_count[index] = np.sum(window)

    return heat_count,heat_label,np.array(divide_rules)

def CreateOverlapMutantData():
    for nsp in nsp_ref:
        print(nsp[0])
        for index,dir in enumerate(folder_path):
            maindata = pd.read_csv("{}{}.csv".format(dir,nsp[0]))
            maindata = maindata.to_numpy()
            ref_seq = nsp[sequence_col_index]
            mutants = np.zeros((maindata.shape[0], 1))
            for i, item in enumerate(maindata):
                founded, m, newm = FindMutant(str(ref_seq).lower(), str(item[sequence_col_index]).lower())
                if (founded == True and np.sum(m)<=15):
                    mutants[i] = np.sum(m)
            try:
                os.mkdir("result/{}".format(nsp[0]))
            except:
                pass

            df = pd.DataFrame({"mutants":mutants[:,0],"ID":maindata[:,1]})
            df.to_csv("result/{}/over.{}.csv".format(nsp[0],sheetname[index]))

def Create_HeatAndPieChart():
    for nsp in nsp_ref:
        print(nsp[0])
        for dir in folder_path:

            csvFile = "{}{}.csv".format(dir,nsp[0])
            if (not os.path.exists(csvFile)):
                continue

            regionName = dir.split('\\')[1].replace("/", "")
            print (regionName)

            dirresult =  "{}/{}/{}".format(prefix_path, nsp[0], regionName)
            if (os.path.exists(dirresult)):
                continue

            maindata = pd.read_csv("{}{}.csv".format(dir,nsp[0]))
            maindata = maindata.to_numpy()
            ref_seq = nsp[sequence_col_index]
            mutants = np.zeros((maindata.shape[0], 1))
            heat_count = np.zeros(10)
            heat_label = None
            overal_find =False

            for i, item in enumerate(maindata):
                founded, m, newm = FindMutant(str(ref_seq).lower(), str(item[sequence_col_index]).lower())
                if (founded == True and np.sum(m)<=15):
                    overal_find=True
                    mutants[i] = np.sum(m)
                    heat_c,heat_label,divide_rules = Find_Heatmap_Data(m)
                    heat_count += heat_c
                if (i==0):
                    heat_c, heat_label, divide_rules = Find_Heatmap_Data(m)

            if (overal_find==True):
                divide_rules = divide_rules*maindata.shape[0]
                heat_count = np.divide(heat_count,divide_rules)

            labels = ["Not Mutant", "One Mutant", "Two Mutant", "Three Mutant", "Four Mutant and more"]

            z0 = np.sum(mutants == 0) / maindata.shape[0]
            z1 = np.sum(mutants == 1) / maindata.shape[0]
            z2 = np.sum(mutants == 2) / maindata.shape[0]
            z3 = np.sum(mutants == 3) / maindata.shape[0]
            z4 = np.sum(mutants >= 4) / maindata.shape[0]

            percetage = np.array([z0, z1, z2, z3, z4])

            dfpie = pd.DataFrame(percetage, index=labels)
            dfheat = pd.DataFrame(heat_count,index=heat_label)

            try:
                os.mkdir("{}/{}".format(prefix_path,nsp[0]))
            except:
                pass
            try:
                regionName = dir.split('\\')[1].replace("/","")
                os.mkdir("{}/{}/{}".format(prefix_path, nsp[0],regionName))
            except:
                pass

            xlsxwriter =  pd.ExcelWriter("{}/{}/{}/pie.xlsx".format(prefix_path,nsp[0],regionName),engine="xlsxwriter")
            dfpie.to_excel(xlsxwriter)
            xlsxwriter.save()

            xlsxwriter =  pd.ExcelWriter("{}/{}/{}/heat.xlsx".format(prefix_path, nsp[0],regionName),engine="xlsxwriter")
            dfheat.to_excel(xlsxwriter)
            xlsxwriter.save()

def CreateTimeLineChart():
    for nsp in nsp_ref:
        print(nsp[0])
        for folder in folder_path:
            csvFile = "{}{}.csv".format(folder,nsp[0])
            fileExit = os.path.exists(csvFile)
            if (os.path.exists(csvFile)==False):
                continue

            maindata = pd.read_csv("{}{}.csv".format(folder, nsp[0]))
            maindata = maindata.to_numpy()

            regionName = folder.split('\\')[1].replace("/", "")
            print (regionName)

            dirresult =  "{}/{}/{}/stacked.xlsx".format(prefix_path, nsp[0], regionName)
            if (os.path.exists(dirresult)):
                continue

            ref_seq = nsp[4]
            date_label =  pd.date_range(start='11/11/2019', periods=30, freq='M')
            distinctdate =  pd.date_range(start='11/11/2019', periods=29, freq='M')
            sequence_count = np.zeros(len(distinctdate) + 1)
            mutants = np.zeros((len(distinctdate)+1, len(ref_seq)))
            total_mutants_new = np.zeros((len(ref_seq),len(aminoacids)))
            for i,item in enumerate(maindata):
                founded, m, newm = FindMutant(str(ref_seq).lower(), str(item[sequence_col_index]).lower())
                total_mutants_new+=newm
                sequencedate = ConvertToDateTime(item[2])

                if (founded == True):
                    if (sequencedate <= distinctdate[0]):
                        mutants[0,  :] += m
                        sequence_count[0] += 1
                    elif (sequencedate > distinctdate[len(distinctdate) - 1]):
                        mutants[len(distinctdate), :] += m
                        sequence_count[len(distinctdate)] += 1
                    else:
                        for dl in range(0, len(distinctdate)):
                            if (sequencedate > distinctdate[dl] and sequencedate <= distinctdate[dl + 1]):
                                mutants[dl + 1, :] += m
                                sequence_count[dl + 1] += 1
                else:
                    if (sequencedate <= distinctdate[0]):
                        sequence_count[0] += 1
                    elif (sequencedate > distinctdate[len(distinctdate) - 1]):
                        sequence_count[len(distinctdate)] += 1
                    else:
                        for dl in range(0, len(distinctdate)):
                            if (sequencedate > distinctdate[dl] and sequencedate <= distinctdate[dl + 1]):
                                sequence_count[dl + 1] += 1
            aa_name = []
            for i,item in enumerate(ref_seq):
                aa_name.append(item.upper()+"({})".format(i+1))
            aa_name =np.array(aa_name)

            #take log2 for all value of dictionary
            final = mutants
            total = np.sum(mutants,axis =0)/maindata.shape[0]
            del mutants
            sortedindex = total.argsort()[::-1][:len(total)]

            array = []

            for i,itemarray in enumerate(final):
                final1 = np.divide(np.transpose(final[i, :]), sequence_count[i])
                final1[np.isnan(final1)] = 0
                final1 = np.round(final1, 6)
                array.append(final1[sortedindex])

            timeline_df = pd.DataFrame(np.transpose(np.array(array))[0:], index=(aa_name[sortedindex])[0:],columns= date_label )

            #make bar plot
            total_mutants_new = np.divide(total_mutants_new,maindata.shape[0])
            total_bar = np.sum(total_mutants_new,axis=1)
            totalsort_index = total_bar.argsort()[::-1][:len(total_bar)]

            new_mutants = total_mutants_new[totalsort_index,:][0:]
            total_bar = np.sum(new_mutants,axis=0)
            col_index = np.array(np.where(total_bar>0.00001))[0]
            if (len(col_index)>0):
                new_mutants = new_mutants[:, col_index]
                barplot_df = pd.DataFrame(new_mutants,index = aa_name[totalsort_index][0:],columns= aminoacids[col_index])
            else:
                new_mutants = new_mutants[:, :]
                barplot_df = pd.DataFrame(new_mutants,index = aa_name[:][0:],columns= aminoacids[:])

            try:
                os.mkdir("{}/{}".format(prefix_path,nsp[0]))
            except:
                pass
            try:
                regionName = folder.split('\\')[1].replace("/","")
                os.mkdir("{}/{}/{}".format(prefix_path, nsp[0],regionName))
            except:
                pass

            xlsxwriter =  pd.ExcelWriter("{}/{}/{}/timeline.xlsx".format(prefix_path,nsp[0],regionName),engine="xlsxwriter")
            timeline_df.to_excel(xlsxwriter)
            xlsxwriter.save()

            xlsxwriter =  pd.ExcelWriter("{}/{}/{}/stacked.xlsx".format(prefix_path, nsp[0],regionName),engine="xlsxwriter")
            barplot_df.to_excel(xlsxwriter)
            xlsxwriter.save()

def Create_DistinctMutant():
    for nsp in nsp_ref:
        print(nsp[0])
        df_mutant = []
        labels = ["ID","Date","Country","Continet","Type"]

        for folder in folder_path:
            csvFile = "{}{}.csv".format(folder,nsp[0])
            if (not os.path.exists(csvFile)):
                continue
            maindata = pd.read_csv("{}{}.csv".format(folder, nsp[0]))
            maindata = maindata.to_numpy()
            ref_seq = nsp[4]
            mutant_type=[]
            mutant_temp = []

            for i, item in enumerate(maindata):
                founded, newm = DistictMutant(str(ref_seq).lower(), str(item[sequence_col_index]).lower())
                if founded == True:
                    temp = [item[1],item[2],item[5],item[6],newm]
                    mutant_type.append(temp)
            df1 = pd.DataFrame(mutant_type,columns= labels)
            df_mutant.append(df1)

            try:
                os.mkdir("{}/{}".format(prefix_path, nsp[0]))
            except:
                pass
            try:
                regionName = folder.split('\\')[1].replace("/", "")
                os.mkdir("{}/{}/{}".format(prefix_path, nsp[0], regionName))
            except:
                pass

            df = pd.DataFrame(mutant_type)
            df.to_csv("{}/{}/{}/mutant.csv".format(prefix_path, nsp[0], regionName))

def FindFirstMutant():
    sheetname = ["Africa","World Wide","Asia","Europe","North America","Oceania","South America"]
    sheetname = ["United States","England","Sudan","Netherlands","France","Russia","Canada","India","Guangdong","Germany","Wuhan","Ecuador",'Kazakhstan','Senegal','Greece','Italy','Belgium','Mexico','Hong Kong','Gambia','Oman','Niger','Sweden','Norway','Japan','Saudi Arabia','New Zealand','Brazil','Austria','Venezuela','Wales','Thailand','Foshan','Denmark']
    folder_path = ["United States","England","Sudan","Netherlands","France","Russia","Canada","India","Guangdong","Germany","Wuhan","Ecuador",'Kazakhstan','Senegal','Greece','Italy','Belgium','Mexico','Hong Kong','Gambia','Oman','Niger','Sweden','Norway','Japan','Saudi Arabia','New Zealand','Brazil','Austria','Venezuela','Wales','Thailand','Foshan','Denmark']
    #reject_nsp =[ "nsp1","nsp2","nsp3","nsp4","nsp5","nsp6","nsp7","nsp8","nsp9","nsp10","nsp11","nsp12","nsp13","nsp14","nsp15","nsp16","spike","envelope","membrane"]
    for nsp in nsp_ref:
        print(nsp[0])
#        if nsp[0] in reject_nsp:
#            continue
        df_firstsaw = []
        for folder,sname in zip(folder_path,sheetname):
            csv_path = "{}{}.csv".format("clean/region/"+ folder+"/", nsp[0])
            if (not os.path.exists(csv_path)):
                continue

            maindata = pd.read_csv("{}{}.csv".format("clean/region/"+ folder+"/", nsp[0]))
            maindata = maindata.to_numpy()
            mydate = maindata[:,2]

            for i in range(len(mydate)):
                mydate[i] = ConvertToDateTime(mydate[i])

            sortedindex = mydate.argsort()
            maindata = maindata[sortedindex,:]
            if (nsp[0] == 'nucleoprotein' ):
                mutants = pd.read_excel(prefix_path+"/{}/{}/stacked.xlsx".format('neoclocapsid', folder))
                mutantscol = mutants.columns[1:]
                mutants = mutants.to_numpy()
            else:
                mutants = pd.read_excel(prefix_path+"/{}/{}/stacked.xlsx".format(nsp[0], folder))
                mutantscol = mutants.columns[1:]
                mutants = mutants.to_numpy()

            firstsaw = np.empty((mutants.shape[0], 5),dtype=np.object)

            for index,mut in enumerate(mutants):
                values = mut[1:]
                maxindex = np.argmax(values)
                mutantrate = np.max(values)
                from1 = mut[0][0].lower()
                position = int(mut[0][2:len(mut[0]) - 1])
                to1 = mutantscol[maxindex].lower()
                founded = False

                for item in maindata:
                    sequence = item[4]
                    country = item[5]
                    continent = item[6]
                    date1 = item[2]
                    at = sequence[position - 1]
                    if (at == to1):
                        founded=True
                        break
                if(founded==True):
                    firstsaw[index,0] = from1+"(" +str(position) +")->"+str(to1)
                    firstsaw[index,1] = mutantrate
                    firstsaw[index,2] = country
                    firstsaw[index,3] = date1
                    firstsaw[index,4] = continent

            df = pd.DataFrame(firstsaw,columns = ["Mutant","Rate","Country","Continent","Date"])
            df_firstsaw.append(df)
            try:
                os.mkdir("result/firstmutant_region/{}".format(nsp[0]))
            except:
                pass

            xlsxwriter =  pd.ExcelWriter("result/firstmutant_region/{}/firstsaw.{}.xlsx".format(nsp[0],folder),engine="xlsxwriter")
            for item,sname in zip(df_firstsaw,sheetname):
                item.to_excel(xlsxwriter,sheet_name=sname)
            xlsxwriter.save()

FindFirstMutant()
#Create_HeatAndPieChart()
#CreateOverlapMutantData()

#CreateTimeLineChart()

#Create_DistinctMutant()
