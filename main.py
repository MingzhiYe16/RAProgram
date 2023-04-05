import csv
import random
import pandas as pd
import numpy
from scipy.special import rel_entr
from scipy.spatial import distance
from scipy.stats import norm
# Global parameter

AmiNoAcids=['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
ProbOfEachClonotype=None
SamplePool=None
realDF=None
joinedDF=None


# IndexnAnBMap=[(0, 0), (2, 1), (1, 2), (1, 1), (0, 2), (1, 0), (2, 2), (0, 1), (3, 1), (1, 3), (2, 3), (3, 2), (2, 0), (3, 3)]
# nAnBWeight=[149, 134, 154, 409, 15, 1, 57, 53, 7, 10, 5, 4, 1, 1]
# IndexnAnBMap=[(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2), (2, 0), (2, 1), (2, 2)]
# nAnBWeight=[0, 12382, 0, 2813, 71977, 2194, 0, 7129, 566]
# IndexnAnBMap=[(1, 1), (2, 2), (4, 3), (2, 1), (3, 2), (1, 2), (0, 1), (1, 0), (0, 2), (4, 2), (2, 3), (3, 3), (3, 4), (4, 4), (3, 1), (5, 3), (1, 3), (2, 4), (5, 5), (5, 4), (4, 5), (0, 3), (2, 0), (3, 5)]
# nAnBWeight=[86223, 8484, 122, 9789, 1301, 5726, 14979, 2342, 188, 78, 692, 654, 71, 37, 33, 10, 114, 33, 3, 9, 3, 1, 1, 1]
# RCC+melanoma w/ new Cell Ranger
# IndexnAnBMap=[(0, 1), (1, 1), (1, 0), (1, 2), (2, 2), (2, 1)]
# nAnBWeight=[30479, 176070, 5612, 6617, 3353, 18511]
# RCC w/ old Cell Ranger:
IndexnAnBMap=[(0, 1), (1, 1), (1, 0), (0, 2), (1, 2), (2, 2), (2, 1), (2, 0), (0, 3), (1, 3), (2, 3), (3, 3), (3, 2), (3, 1), (3, 0), (0, 4), (1, 4), (2, 4), (3, 4), (4, 4), (4, 3), (4, 2), (4, 1), (4, 0), (1, 5), (2, 5), (3, 5), (4, 5), (5, 5), (5, 4), (5, 3), (5, 2), (5, 1), (1, 6), (2, 6), (3, 6), (5, 6), (6, 5), (6, 4), (6, 2), (7, 6), (2, 7), (3, 7)]
nAnBWeight=[9469, 50038, 3577, 2245, 25671, 10315, 13407, 102, 90, 1504, 1538, 629, 1171, 516, 5, 5, 74, 146, 80, 34, 105, 126, 30, 1, 11, 12, 11, 9, 3, 7, 14, 10, 1, 1, 5, 5, 2, 1, 2, 1, 1, 1, 1]


def generateSamplePool(PoolSize):
    global SamplePool,realDF,joinedDF
    joinedDF = None
    # Generate sample pool
    def getCDR3():
        res=random.choices(AmiNoAcids,k=15)
        return ''.join(res)
    Original_Alpha_Pool=[getCDR3() for i in range(int(PoolSize/2))]
    Original_Beta_Pool=[getCDR3() for i in range(int(PoolSize/2))]
    SamplePoolAlpha=random.choices(Original_Alpha_Pool,k=PoolSize)
    SamplePoolBeta=random.choices(Original_Beta_Pool,k=PoolSize)
    SamplePool=list(zip(SamplePoolAlpha,SamplePoolBeta))


    return

def getSample(PoolSize,SampleSize,BaseNumber,errorProb=0.01):
    global ProbOfEachClonotype, realDF
    # ProbOfEachClonotype=[1 for i in range(PoolSize)]

    if BaseNumber>0:
        ProbOfEachClonotype = [random.uniform(0.5,0.95) ** i for i in range(PoolSize)]
    else:
        ProbOfEachClonotype = [BaseNumber ** i + 0.999 ** i for i in range(PoolSize)]
    ProbOfEachClonotype = [x if x > 0.1 ** 6 else 0.1 ** 6 for x in ProbOfEachClonotype]

    ProbSum = sum(ProbOfEachClonotype)
    ProbOfEachClonotype = [i / ProbSum for i in ProbOfEachClonotype]

    real_clonotypes = [x[0] + ' ' + x[1] for x in SamplePool]
    real_abundance = ProbOfEachClonotype[::]
    names = ("clonotypes", "abundance")
    realDF = pd.DataFrame(list(zip(real_clonotypes, real_abundance)),
                          columns=names)
    realDF = realDF.groupby(['clonotypes']).sum()
    realDF = realDF.reset_index()
    NumOfChains=random.choices(IndexnAnBMap, weights=nAnBWeight,k=SampleSize)

    def getChainTypes(i):
        res=[]
        for j in range(NumOfChains[i][0]):
            res.append('TRA')
        for j in range(NumOfChains[i][1]):
            res.append('TRB')
        return res

    chains=[getChainTypes(i) for i in range(SampleSize)]
    chains = [item for sublist in chains for item in sublist]

    def getBarcode(i):
        res=[]
        for j in range(sum(NumOfChains[i])):
            res.append('cell_'+str(i+1))

        return res

    barcodes=[getBarcode(i) for i in range(SampleSize)]
    barcodes = [item for sublist in barcodes for item in sublist]

    def getCDR3sForCertainClonotype(i):
        result=[[],[]]
        numberA,numberB=NumOfChains[i][0],NumOfChains[i][1]
        indexInPool = numpy.random.choice(range(PoolSize),size=3,replace=False, p=ProbOfEachClonotype)
        usedA,usedB=set(),set()
        if numberA >= 1:
            if random.uniform(0,1) >= errorProb:
                result[0].append(SamplePool[indexInPool[0]][0])
                usedA.add(indexInPool[0])
            else:
                result[0].append(SamplePool[indexInPool[1]][0])
                usedA.add(indexInPool[1])
        numberA-=1

        if numberB >= 1:
            if random.uniform(0,1) >= errorProb:
                result[1].append(SamplePool[indexInPool[0]][1])
                usedB.add(indexInPool[0])
            else:
                result[1].append(SamplePool[indexInPool[2]][1])
                usedB.add(indexInPool[2])
        numberB-=1

        while numberA > 0:
            indexInPool=random.choices(range(PoolSize),k=1,weights=ProbOfEachClonotype)
            indexInPool=indexInPool[0]
            while indexInPool in usedA:
                indexInPool = random.choices(range(PoolSize),k=1,weights=ProbOfEachClonotype)
                indexInPool = indexInPool[0]
            usedA.add(indexInPool)
            result[0].append(SamplePool[indexInPool][0])
            numberA-=1

        while numberB > 0:
            indexInPool=random.choices(range(PoolSize),k=1,weights=ProbOfEachClonotype)
            indexInPool = indexInPool[0]
            while indexInPool in usedB:
                indexInPool = random.choices(range(PoolSize),k=1,weights=ProbOfEachClonotype)
                indexInPool = indexInPool[0]
            usedB.add(indexInPool)
            result[1].append(SamplePool[indexInPool][1])
            numberB-=1

        return result[0]+result[1]

    CDR3=[getCDR3sForCertainClonotype(i) for i in range(SampleSize)]
    CDR3 = [item for sublist in CDR3 for item in sublist]

    num_of_chain=len(chains)
    IS_CELL=['TRUE']*num_of_chain
    contig_id = ['CHAR'] * num_of_chain
    high_confidence = ['TRUE'] * num_of_chain
    LENGTH = [697] * num_of_chain
    v_gene = ['CHAR'] * num_of_chain
    d_gene = ['CHAR'] * num_of_chain
    c_gene = ['CHAR'] * num_of_chain
    j_gene = ['CHAR'] * num_of_chain
    full_length = ['TRUE'] * num_of_chain
    productive = ['TRUE'] * num_of_chain
    cdr3_nt = ['CHAR'] * num_of_chain
    reads = [44152] * num_of_chain
    umis = [15] * num_of_chain
    raw_clonotype_id = ['CHAR'] * num_of_chain
    TYPE0 = ['TCR'] * num_of_chain
    clonotype_tcr = ['CHAR'] * num_of_chain
    raw_consensus_id = ['CHAR'] * num_of_chain
    SAMPLE = ['CHAR'] * num_of_chain


    sample=list(zip(barcodes,IS_CELL,contig_id,high_confidence,LENGTH,chains,v_gene,d_gene,j_gene,c_gene \
             ,full_length,productive,CDR3,cdr3_nt, reads, umis, raw_clonotype_id, raw_consensus_id, \
             TYPE0, clonotype_tcr,SAMPLE))

    names= ( "barcode","is_cell","contig_id","high_confidence","length","chain","v_gene","d_gene","j_gene","c_gene" \
             ,"full_length","productive","cdr3", "cdr3_nt", "reads", "umis", "raw_clonotype_id", "raw_consensus_id", \
             "type", "clonotype_tcr","sample")
    outputDF = pd.DataFrame(sample,
                      columns=names)

    return outputDF



def updateMergedDataframe(clonotypes,abundance):
    global joinedDF
    names=("clonotypes","abundance")
    stimulated = pd.DataFrame(list(zip(clonotypes,abundance)),
                            columns=names)
    sum_abundance=sum(stimulated['abundance'])
    stimulated['abundance']=stimulated['abundance']/sum_abundance
    if joinedDF is None:
        joinedDF=realDF.copy()
    joinedDF=joinedDF.join(stimulated.set_index(['clonotypes']),how='outer', on=['clonotypes'],lsuffix='1', rsuffix='2')
    joinedDF = joinedDF.fillna(0)
    return

def getIntermediateResult(Times):
    global joinedDF
    names=[['EM'+str(i+1),'UN'+str(i+1)] for i in range(Times)]
    names=[a for b in names for a in b]
    names = ['clonotypes', 'real'] + names
    joinedDF.set_axis(names, axis=1, inplace=True)


    return joinedDF

def clearStorage():
    global ProbOfEachClonotype,SamplePool,realDF,joinedDF
    ProbOfEachClonotype = None
    SamplePool = None
    realDF = None
    joinedDF = None

def getDistance(clonotypes,abundance):
    names=("clonotypes","abundance")
    stimulated = pd.DataFrame(list(zip(clonotypes,abundance)),
                            columns=names)
    sum_abundance=sum(stimulated['abundance'])
    stimulated['abundance']=stimulated['abundance']/sum_abundance
    joinedDF=realDF.join(stimulated.set_index(['clonotypes']),how='outer', on=['clonotypes'],lsuffix='1', rsuffix='2')
    joinedDF = joinedDF.fillna(0)
    # return joinedDF
    return sum(abs(joinedDF['abundance1']-joinedDF['abundance2']))

def getDistanceTwoSamples(clonotypes1,abundance1,clonotypes2,abundance2):
    names=("clonotypes","abundance")
    stimulated1 = pd.DataFrame(list(zip(clonotypes1,abundance1)),
                            columns=names)
    stimulated2 = pd.DataFrame(list(zip(clonotypes2,abundance2)),
                            columns=names)
    sum_abundance1=sum(stimulated1['abundance'])
    stimulated1['abundance']=stimulated1['abundance']/sum_abundance1
    sum_abundance2=sum(stimulated2['abundance'])
    stimulated2['abundance']=stimulated2['abundance']/sum_abundance2

    joinedDF=stimulated1.join(stimulated2.set_index(['clonotypes']),how='outer', on=['clonotypes'],lsuffix='1', rsuffix='2')
    joinedDF = joinedDF.fillna(0)
    # return joinedDF
    return sum(abs(joinedDF['abundance1']-joinedDF['abundance2']))

def getKLDivergence(clonotypes,abundance):
    names=("clonotypes","abundance")
    stimulated = pd.DataFrame(list(zip(clonotypes,abundance)),
                            columns=names)
    sum_abundance=sum(stimulated['abundance'])
    stimulated['abundance']=stimulated['abundance']/sum_abundance
    joinedDF=realDF.join(stimulated.set_index(['clonotypes']),how='outer', on=['clonotypes'],lsuffix='1', rsuffix='2')
    joinedDF = joinedDF.fillna(1e-20)
    # return joinedDF
    return sum(rel_entr(joinedDF['abundance2'],joinedDF['abundance1']))

def getJSDivergence(clonotypes,abundance):
    names=("clonotypes","abundance")
    stimulated = pd.DataFrame(list(zip(clonotypes,abundance)),
                            columns=names)
    sum_abundance=sum(stimulated['abundance'])
    stimulated['abundance']=stimulated['abundance']/sum_abundance
    joinedDF=realDF.join(stimulated.set_index(['clonotypes']),how='outer', on=['clonotypes'],lsuffix='1', rsuffix='2')
    joinedDF = joinedDF.fillna(0)
    # return joinedDF
    return distance.jensenshannon(joinedDF['abundance2'],joinedDF['abundance1'])
########################################################################
#############       Find rules for naming clonotypes       #############
########################################################################

# import pandas as pd
# df = pd.read_csv ('C:/Users/yemin/Downloads/raw_data/raw_data\S1_T/filtered_contig_annotations.csv')
# df1=df[["raw_clonotype_id","chain","cdr3","exact_subclonotype_id","barcode"]]
#
# df2=df1[df1["raw_clonotype_id"]=='clonotype3606']
# df2.head(50)
#
# df3=df1[df1["cdr3"]=='CAFMKHTYPQGGSEKLVF']
# df3.head(50)
#
# s1=df.groupby('barcode')['raw_clonotype_id'].apply(list)
# s2=df.groupby('barcode')['cdr3'].apply(list)
# df4=pd.DataFrame()
# df4['raw_clonotype_id']=s1
# df4['cdr3']=s2
# df4['length']=1
# def filterunique(x):
#     x[2]=len(x[1])
#     return x
# df4=df4.apply(lambda x:filterunique(x),axis=1)
# df4=df4[df4.length==1]
# def delist(x):
#     x[1]=x[1][0]
#     x[0]=x[0][0]
#     return x
# df4=df4.apply(lambda x:delist(x),axis=1)
# s3=df4.groupby('cdr3')['raw_clonotype_id'].apply(list)
# df5=pd.DataFrame()
# df5["clonotype"]=s3
# df5['length']=1
# def prep1(x):
#     x[0]=set(x[0])
#     x[1]=len(x[0])
#     return x
# df5=df5.apply(lambda x:prep1(x),axis=1)
# df5=df5[df5.length!=1]

#
# df1.sort_values(by=['raw_clonotype_id',"column_index"])
# df1.head(40)



###################################################################################
#############      Explore distribution of number of TRA and TRB      #############
###################################################################################

# import pandas as pd
# import glob
# import os
# from collections import Counter
# path = r'/Users/mingzhiye/RProjects/RAProgram/Melanoma_JEM_TCR_DataPrep' # use your path
# all_files = glob.glob(os.path.join(path , "*/*.csv"))
#
# li = []
#
# for filename in all_files:
#     df = pd.read_csv(filename, index_col=None, header=0)
#     li.append(df)
#
# df = pd.concat(li, axis=0, ignore_index=True)
# s1=df.groupby('barcode')['chain'].apply(list)
# df2=pd.DataFrame()
# df2['chain']=s1
# df2['distribution']=1
# def getDistribution(x):
#     distr=[0,0]
#     for s in x[0]:
#         if s=='TRA':
#             distr[0]+=1
#         elif s=='TRB':
#             distr[1]+=1
#     x[1]=distr
#     return x
# df2=df2.apply(lambda x:getDistribution(x),axis=1)
# distr=list(df2["distribution"])
# distr=[tuple(x) for x in distr]
# c=Counter(distr)
