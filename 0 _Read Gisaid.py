from Bio.Seq import Seq

import pandas as pd
import numpy as np

class VirusList():
    def __init__(self):
       self.proteinID=""
       self.proteinSequence =  ""
class Protein:
    def ReadProteinFasta(self, filename):
        file = open(filename, "r")
        lines = file.readlines()
        index = 0
        print("Lines Count: {}".format(len(lines)))
        virusList =[]
        while (index < len(lines) - 1):
            if (lines[index].find(">") >= 0):
                org = VirusList()
                org.proteinID = lines[index].replace("\n", "").split('|')[0].replace(">", "").split(".")[0]
                sequence = ""
                index += 1
                total = []
                while (lines[index].find(">") < 0 and index < len(lines) - 1):
                    line = lines[index].replace("\n", "")
                    sequence += line
                    index += 1
                    if (len(sequence) >= 500):
                        total.append(sequence.lower())
                        sequence = ""
                total.append(sequence)
                org.proteinSequence = "".join(total)
                virusList.append(org)
        print("Protein Count = {}".format(len(virusList)))
        return virusList
    def SaveToExel(self,data):
        data = np.array(data)
        newdata = np.empty((data.shape[0], 2), dtype=np.object)
        for i,item in enumerate(data):
            newdata[i,0] = str(item.proteinID).strip().lower()
            newdata[i,1] = str(item.proteinSequence).strip().lower()
        df = pd.DataFrame(newdata,columns=["ProteinID","Sequence"])
        df.to_excel('Data/1_Step (Protein_Read)/maindata.xlsx',index=False)
i=1
filepath = 'data/info/covid_sequence.fasta'

nsp1file  = open('data/nsp1.fasta',"w")
nsp2file  = open('data/nsp2.fasta',"w")
nsp3file  = open('data/nsp3.fasta',"w")
nsp4file  = open('data/nsp4.fasta',"w")
nsp5file  = open('data/nsp5.fasta',"w")
nsp6file  = open('data/nsp6.fasta',"w")
nsp7file  = open('data/nsp7.fasta',"w")
nsp8file  = open('data/nsp8.fasta',"w")
nsp9file  = open('data/nsp9.fasta',"w")
nsp10file = open('data/nsp10.fasta',"w")
nsp11file = open('data/nsp11.fasta',"w")
nsp12file = open('data/nsp12.fasta',"w")
nsp13file = open('data/nsp13.fasta',"w")
nsp14file = open('data/nsp14.fasta',"w")
nsp15file = open('data/nsp15.fasta',"w")
nsp16file = open('data/nsp16.fasta',"w")

spikefile = open('data/spike.fasta',"w")
Efile = open('data/envelope.fasta',"w")
Mfile = open('data/membrane.fasta',"w")
Nfile = open('data/nucleoprotein.fasta',"w")

orf3afile = open('data/orf3a.fasta',"w")
orf6file =  open('data/orf6.fasta',"w")
orf7afile = open('data/orf7a.fasta',"w")
orf7bfile = open('data/orf7b.fasta',"w")
orf8file =  open('data/orf8.fasta',"w")
orf10file = open('data/orf10.fasta',"w")
orf9cfile = open('data/orf9c.fasta',"w")
orf9bfile = open('data/orf9b.fasta',"w")


nsp1_count =0
nsp2_count =0
nsp3_count =0
nsp4_count =0
nsp5_count =0
nsp6_count =0
nsp7_count =0
nsp8_count =0
nsp9_count =0
nsp10_count =0
nsp11_count =0
nsp12_count =0
nsp13_count =0
nsp14_count =0
nsp15_count =0
nsp16_count =0

spike_count =0
efile_count =0
mfile_count =0
nfile_count =0

orf3a_count = 0
orf6_count = 0
orf7a_count = 0
orf7b_count = 0
orf8_count = 0
orf9b_count = 0
orf9c_count = 0
orf10_count = 0

i=1
with open(filepath) as fp:
    while (1==1):
        i+=1
        print(i)

        line1 = fp.readline()
        if (not line1):
            break

        line2 = fp.readline()
        header = line1.split("|")

        if  (header[0] == ">NSP1"):
            nsp1file.write(line1)
            nsp1file.write(str(line2).lower())
            nsp1_count+=1
        if  (header[0] == ">NSP2"):
            nsp2file.write(line1)
            nsp2file.write(str(line2).lower())
            nsp2_count+=1
        if  (header[0] == ">NSP3"):
            nsp3file.write(line1)
            nsp3file.write(str(line2).lower())
            nsp3_count+=1
        if  (header[0] == ">NSP4"):
            nsp4file.write(line1)
            nsp4file.write(str(line2).lower())
            nsp4_count+=1
        if  (header[0] == ">NSP5"):
            nsp5file.write(line1)
            nsp5file.write(str(line2).lower())
            nsp5_count+=1
        if  (header[0] == ">NSP6"):
            nsp6file.write(line1)
            nsp6file.write(str(line2).lower())
            nsp6_count+=1
        if  (header[0] == ">NSP7"):
            nsp7file.write(line1)
            nsp7file.write(str(line2).lower())
            nsp7_count+=1
        elif (header[0] == ">NSP8"):
            nsp8file.write(line1)
            nsp8file.write(line2)
            nsp8_count += 1
        elif (header[0] == ">NSP9"):
            nsp9file.write(line1)
            nsp9file.write(line2)
            nsp9_count += 1
        elif (header[0] == ">NSP10"):
            nsp10file.write(line1)
            nsp10file.write(line2)
            nsp10_count += 1
        elif (header[0] == ">NSP11"):
            nsp11file.write(line1)
            nsp11file.write(line2)
            nsp11_count += 1
        elif (header[0] == ">NSP12"):
            nsp12file.write(line1)
            nsp12file.write(line2)
            nsp12_count += 1
        elif (header[0] == ">NSP13"):
            nsp13file.write(line1)
            nsp13file.write(line2)
            nsp13_count += 1
        elif (header[0] == ">NSP14"):
            nsp14file.write(line1)
            nsp14file.write(line2)
            nsp14_count += 1
        elif (header[0] == ">NSP15"):
            nsp15file.write(line1)
            nsp15file.write(line2)
            nsp15_count += 1
        elif (header[0] == ">NSP16"):
            nsp16file.write(line1)
            nsp16file.write(line2)
            nsp16_count += 1
        elif (header[0] == ">Spike"):
            spikefile.write(line1)
            spikefile.write(line2)
            spike_count += 1
        elif (header[0] == ">E"):
            Efile.write(line1)
            Efile.write(line2)
            efile_count += 1
        elif (header[0] == ">N"):
            Nfile.write(line1)
            Nfile.write(line2)
            nfile_count += 1
        elif (header[0] == ">M"):
            Mfile.write(line1)
            Mfile.write(line2)
            mfile_count += 1
        elif (header[0] == ">NS3"):
            orf3afile.write(line1)
            orf3afile.write(line2)
            orf3a_count += 1
        elif (header[0] == ">NS6"):
            orf6file.write(line1)
            orf6file.write(line2)
            orf6_count += 1
        elif (header[0] == ">NS7a"):
            orf7afile.write(line1)
            orf7afile.write(line2)
            orf7a_count += 1
        elif (header[0] == ">NS7b"):
            orf7bfile.write(line1)
            orf7bfile.write(line2)
            orf7b_count += 1
        elif (header[0] == ">NS8"):
            orf8file.write(line1)
            orf8file.write(line2)
            orf8_count += 1
        elif (header[0] == ">NS10"):
            orf10file.write(line1)
            orf10file.write(line2)
            orf10_count += 1
        elif (header[0] == ">NS9b"):
            orf9bfile.write(line1)
            orf9bfile.write(line2)
            orf9b_count += 1
        elif (header[0] == ">NS9c"):
            orf9cfile.write(line1)
            orf9cfile.write(line2)
            orf9c_count += 1
nsp1file.close()
nsp2file.close()
nsp3file.close()
nsp4file.close()
nsp5file.close()
nsp6file.close()
nsp7file.close()
nsp8file.close()
nsp9file.close()
nsp10file.close()
nsp11file.close()
nsp12file.close()
nsp13file.close()
nsp14file.close()
nsp15file.close()
nsp16file.close()

spikefile.close()
Efile.close()
Nfile.close()
Mfile.close()

orf3afile.close()
orf6file.close()
orf7afile.close()
orf7bfile.close()
orf8file.close()
orf9bfile.close()
orf9cfile.close()
orf10file.close()



print ("Nsp 1 :{} ".format(nsp1_count))
print ("Nsp 2 :{} ".format(nsp2_count))
print ("Nsp 3 :{} ".format(nsp3_count))
print ("Nsp 4 :{} ".format(nsp4_count))
print ("Nsp 5 :{} ".format(nsp5_count))
print ("Nsp 6 :{} ".format(nsp6_count))
print ("Nsp 7 :{} ".format(nsp7_count))
print ("Nsp 8 :{} ".format(nsp8_count))
print ("Nsp 9 :{} ".format(nsp9_count))
print ("Nsp 10 :{} ".format(nsp10_count))
print ("Nsp 11 :{} ".format(nsp11_count))
print ("Nsp 12 :{} ".format(nsp12_count))
print ("Nsp 13 :{} ".format(nsp13_count))
print ("Nsp 14 :{} ".format(nsp14_count))
print ("Nsp 15 :{} ".format(nsp15_count))
print ("Nsp 16 :{} ".format(nsp16_count))

print ("spike :{} ".format(spike_count))
print ("Nucleocapcid :{} ".format(nfile_count))
print ("Envelope :{} ".format(efile_count))
print ("Membrane :{} ".format(mfile_count))

print ("ORF 3a :{} ".format(orf3a_count))
print ("ORF 6 :{} ".format(orf6_count))
print ("ORF 7a :{} ".format(orf7a_count))
print ("ORF 7a :{} ".format(orf7b_count))
print ("ORF 8 :{} ".format(orf8_count))
print ("ORF 9b :{} ".format(orf9b_count))
print ("ORF 9c :{} ".format(orf9c_count))
print ("ORF 10 :{} ".format(orf10_count))

nspname = ["Nsp1","Nsp2","Nsp3","Nsp4","Nsp5","Nsp6","Nsp7","Nsp8","Nsp9","Nsp10","Nsp11","Nsp12","Nsp13","Nsp14","Nsp15","nsp16","Spike","Envelope","Nucleprotein","Membrane","Orf3a","Orf6","Orf7a","Orf7b","Orf8","Orf9b","Orf9c","Orf10"]
gene_count = [nsp1_count,nsp2_count,nsp3_count,nsp4_count,nsp5_count,nsp6_count,nsp7_count,nsp8_count,nsp9_count,nsp10_count,nsp11_count,nsp12_count,nsp13_count,nsp14_count,nsp15_count,nsp16_count,spike_count,nfile_count,efile_count,mfile_count,orf3a_count,orf6_count,orf7a_count,orf7b_count,orf8_count,orf9b_count,orf9c_count,orf10_count]

df = pd.DataFrame({"NspName":nspname,"Count":gene_count})
df.to_excel("data/Statistic.xlsx")
