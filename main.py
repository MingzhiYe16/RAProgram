import csv
import random
import pandas as pd

# Global parameter
PoolSize=2000
SampleSize=5000
errorProb=0.01
AmiNoAcids=['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']

# Generate sample pool
def getCDR3():
    res=random.choices(AmiNoAcids,k=15)
    return ''.join(res)
Original_Alpha_Pool=[getCDR3() for i in range(PoolSize)]
Original_Beta_Pool=[getCDR3() for i in range(PoolSize)]
SamplePoolAlpha=random.choices(Original_Alpha_Pool,k=PoolSize)
SamplePoolBeta=random.choices(Original_Beta_Pool,k=PoolSize)
SamplePool=list(zip(SamplePoolAlpha,SamplePoolBeta))

IndexnAnBMap=[]
with open('IndexnAnBMap.csv', 'r') as csvfile:
    reader = csv.reader(csvfile, skipinitialspace=True)
    for row in reader:
        IndexnAnBMap.append(row)
IndexnAnBMap=[[int(x[0]),int(x[1])] for x in IndexnAnBMap]


real_clonotypes=[x[0]+' '+x[1] for x in SamplePool]
real_abundance=[1/PoolSize for i in range(PoolSize)]
names=("clonotypes","abundance")
realDF=pd.DataFrame(list(zip(real_clonotypes,real_abundance)),
                            columns=names)
realDF=realDF.groupby(['clonotypes']).sum()
realDF=realDF.reset_index()
def getSample():
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
        indexInPool = random.sample(range(PoolSize), 3)
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
            indexInPool=random.randint(0, PoolSize-1)
            while indexInPool in usedA:
                indexInPool = random.randint(0, PoolSize - 1)
            usedA.add(indexInPool)
            result[0].append(SamplePool[indexInPool][0])
            numberA-=1

        while numberB > 0:
            indexInPool=random.randint(0, PoolSize-1)
            while indexInPool in usedB:
                indexInPool = random.randint(0, PoolSize - 1)
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



