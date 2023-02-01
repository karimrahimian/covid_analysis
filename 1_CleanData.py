import pandas as pd
import numpy as np
import os

distinct_country = dict()

def findContinetName(country,continent_ref):
    for i, item in enumerate(continent_ref):
        if (str(item[1]).lower() == str(country).lower()):
            return item[2].lower()
    return "NA"
def correctCountry(un_correct_country):
    un_correct_country = un_correct_country.strip()
    for country in correct_country_names:
        names = country.split(",")
        founded = False
        for name in names:
            temp = name.strip()
            temp = temp.replace("#"," ")
            if (temp == un_correct_country):
                founded=True
                break
        if (founded==True):
            correct_name = names[0].replace("#","")
            return correct_name
    return un_correct_country
def cleannsp(nspfilename,length):
    i = 1
    notfound = 0
    data = []

    mylength =dict()
    sequence_with_x = 0
    sequence_not_human = 0
    notlengthmatch = 0
    with open(nspfilename) as fp:
        while (1==1):
            #print(i)
            line1 = fp.readline()
            if (not line1):
                break
            head =line1.split("|")
            line2 = fp.readline().replace("\n" ,"").replace("*","")
            #print (len(line2))
            try:
                mylength[len(line2)] += 1
            except:
                mylength[len(line2)] = 1

            i += 1
            if (len(line2)!=length):
                notfound+=1
                notlengthmatch+=1
                continue

            if (str(line2).lower().find("x")>=0):
                notfound += 1
                sequence_with_x+=1
                continue

            if (head[6]!= 'Human'):
                notfound += 1
                sequence_not_human+=1
                continue

            count = head[1].split("/")
            serial1 = ">"+count[2]
            serial2 = head[3]

            country = str(count[1]).strip().replace("_"," ")
            country = correctCountry(country)

            genotype = head[6]
            continent_name = findContinetName(country,continent_ref)

            row = [ serial1 , serial2,  head[2].replace("\n" ,""),  length, line2.lower().replace("\n" ,"").replace("*","") , country , continent_name ,genotype]
            try:
                distinct_country[country]=1
            except Exception as ex:
                print(ex)

            data.append(row)

    report_file.write("********************************************************************************************************\n")
    report_file.write("{} , Not Found : {} , Total : {} , Founded : {} , Sequence ContainsX: {} , Not Human : {}\n" .format(nspfilename, notfound,i , i-notfound,sequence_with_x,sequence_not_human))
    for key in mylength:
        report_file.write("{}:{} , ".format(key,mylength[key]))
    report_file.write("\n")
    statistic = [nspfilename,notfound,i,i-notfound,sequence_with_x,sequence_not_human,notlengthmatch]

#    print ("{} , Not Found : {} , Total : {} , Founded : {} , Sequence ContainsX: {} , Not Human : {}, Not Length Match : {}" .format(nspfilename, notfound,i , i-notfound,sequence_with_x,sequence_not_human,notlengthmatch))
#    print (mylength)

    df = pd.DataFrame(data, columns=["ID1","ID2","Date", "Length","Sequence","Country","Continent" , "GenType"])
    filename = os.path.basename(nspfilename)
    df.to_csv('clean/total/{}.csv'.format(filename.replace(".fasta","")), index=False)
    return statistic

continent_ref = pd.read_excel("clean/AllCountry.xlsx", sheet_name='AllCountry')
continent_ref = continent_ref.to_numpy()

ref_seq = pd.read_excel('data/info/ref_sequence.xlsx')
ref = ref_seq.to_numpy()

correct_country_file = open('clean/correct_country.txt', "r")
correct_country_names = correct_country_file.readlines()
correct_country_file.close()

report_file = open("clean/clean_report.txt","w")
statistics = []
for item in ref:
    try:
        print (item[0])
        stat = cleannsp("data/{}.fasta".format(item[0]) ,item[3])
        statistics.append(stat)
    except Exception  as ex:
        print("Error for {}".format(ex))
        continue
report_file.close()

df = pd.DataFrame(np.array(statistics), columns=["NSP", "Not Found", "Total", "Remained", "Sequence contain X", "Not human", "Not Match Length"])
df.to_csv('clean/statistic.csv', index=False)

country_file = open("clean/country.csv","w")
for item in distinct_country.keys():
    country_file.write("{}\n".format(item))
country_file.close()



