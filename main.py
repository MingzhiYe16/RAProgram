import csv
import random
import pandas as pd
import numpy
from scipy.special import rel_entr
from scipy.spatial import distance
# Global parameter

AmiNoAcids=['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
IndexnAnBMap=[]
ProbOfEachClonotype=None
SamplePool=None
realDF=None
joinedDF=None

with open('IndexnAnBMap.csv', 'r') as csvfile:
    reader = csv.reader(csvfile, skipinitialspace=True)
    for row in reader:
        IndexnAnBMap.append(row)
IndexnAnBMap=[[int(x[0]),int(x[1])] for x in IndexnAnBMap]

def generateSamplePool(PoolSize):
    global ProbOfEachClonotype,SamplePool,realDF,joinedDF
    joinedDF = None
    # Generate sample pool
    def getCDR3():
        res=random.choices(AmiNoAcids,k=15)
        return ''.join(res)
    Original_Alpha_Pool=[getCDR3() for i in range(PoolSize)]
    Original_Beta_Pool=[getCDR3() for i in range(PoolSize)]
    SamplePoolAlpha=random.choices(Original_Alpha_Pool,k=PoolSize)
    SamplePoolBeta=random.choices(Original_Beta_Pool,k=PoolSize)
    SamplePool=list(zip(SamplePoolAlpha,SamplePoolBeta))

    ProbOfEachClonotype=[1/(i+1) for i in range(PoolSize)]
    ProbSum=sum(ProbOfEachClonotype)
    ProbOfEachClonotype=[i / ProbSum for i in ProbOfEachClonotype]

    real_clonotypes=[x[0]+' '+x[1] for x in SamplePool]
    real_abundance=ProbOfEachClonotype[::]
    names=("clonotypes","abundance")
    realDF=pd.DataFrame(list(zip(real_clonotypes,real_abundance)),
                                columns=names)
    realDF=realDF.groupby(['clonotypes']).sum()
    realDF=realDF.reset_index()
    return

def getSample(PoolSize,SampleSize,errorProb=0.01):
    NumOfChains=random.choices(IndexnAnBMap,k=SampleSize)

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

def getDistanceUN_EM(clonotypes1,abundance1,clonotypes2,abundance2):
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
