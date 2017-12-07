""" 
#####################################################################################################
This module is for finding genes that are significantly correlated with phenotypes or other genotypes


Inputs:

     
Outputs:


Dan Olson
11-4-2016
#####################################################################################################
"""

#----------------------------------------------------------------------------------------------------
# imports
#----------------------------------------------------------------------------------------------------
from scipy.special import comb
from math import pow
import numpy as np
import pandas as pd

#----------------------------------------------------------------------------------------------------
# constants and globals
#----------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------
# functionso
#----------------------------------------------------------------------------------------------------
def findGenesCorrelatedWithFermData(agFermDict, fermProdList, omAnno, minCnt = 3, sigVal = 0.05, minReadFrac = 0.7):
    """
    find genes with mutations that are significantly correlated with changes in fermentation data
    compare origin mutations with parent-child fermentation data ratios
    agFermDict is a dictionary of parent-child fermentation data from getAllAgFermData (get_ferm_data_v4.py)
    """
    minReadFracFilter = omAnno['readFrac'] > minReadFrac
    dfList = []
    for prod in fermProdList:
        print('\nanalyzing product: ',prod)
        fermDf = agFermDict[prod].copy()
        minLogRatio = findLogRatioFromMinCnt(fermDf, sigVal, minCnt)
        print('change needed for significance = {:.0f}%'.format(100*pow(minLogRatio, 2)))
        hasPhenotype = abs(fermDf['logRatio']) >= minLogRatio
        phenotypeFraction = len(fermDf[hasPhenotype]) / len(fermDf)
        strainWithFerm = fermDf['StrainID'].unique().tolist()
        rightStrain = omAnno['Strain'].isin(strainWithFerm)
        cdsList = omAnno['Locus Tag'].unique().tolist()
    
        resultList = []
        for cds in cdsList:
            #cds = 'Clo1313_1185'
    
            rightLocusTag = omAnno['Locus Tag'] == cds
            strainWithMutList = omAnno.loc[rightStrain & minReadFracFilter & rightLocusTag, 'Strain'].unique().tolist()
            strainWithPhenoList = fermDf.loc[hasPhenotype, 'StrainID'].tolist()
            mutAndPheno = list(set(strainWithMutList).intersection(set(strainWithPhenoList)))
            mutNoPheno = list(set(strainWithMutList).difference(set(strainWithPhenoList)))
    
            sig = sigGene(phenotypeFraction, len(mutAndPheno),len(strainWithMutList) )
            resultList.append([cds, sig, mutAndPheno, mutNoPheno, prod])
            #print('cds={}, significance={:.2f}'.format(cds, sig))
    
        sigCdsDf = pd.DataFrame(resultList, columns = ['CDS', 'Significance', 'MutWithPheno', 'MutNoPheno', 'Phenotype']).sort_values('Significance')
        prodResult = sigCdsDf.loc[sigCdsDf['Significance'] < 2*sigVal, :]
        dfList.append(prodResult)
    
    result = pd.concat(dfList)
    return result


def sigGene(sigFrac, sigCount, totalCount):
    """
    determine the probability that n mutations is a given CDS are significant
    "sigFrac" is the fraction of the population with the phenotype of interest
    "sigCount" is the number of strains with the phenotype of interest that also have a 
         mutation in the CDS of interest
    "totalCount" is the total number of strains with a mutation in the CDS of interest
    The result is the probability of an outcome at least this extreme arising by chance
    """
    probList = [] # probabilities for each scenario
    for i in range(totalCount - sigCount + 1):
        prob = pow(sigFrac,sigCount) * pow(1-sigFrac, totalCount-sigCount) * comb(totalCount, sigCount)
        probList.append(prob)
        #print(prob, ' ', i)
        # increase sigCount to take into account the "or more" part of the algorithm
        sigCount += 1
    
    return sum(probList)
    
    
def sigVsMutCnt(minSignificance, plotResult = False):
    """ 
    determine the minimum number of mutations needed to give a significant value
    for a given significance fraction.  For example, if the minSignificance is set to 0.05,
    and you wanted to be able to detect significant results with just 2 mutations 
    (i.e. 2 strains with the phenotype have the mutation and no strains without the 
    phenotype have the mutation), then you would need to choose a phenotype cutoff
    such that the sigFrac is 0.22 or lower
    """
    #minSignificance = 0.05
    stepSize = 0.001
    resultList = []
    for mutCnt in range(1,8):
        for sigFrac in np.arange(0,1,stepSize):
            sig = sigGene(sigFrac, mutCnt, mutCnt)
            if sig > minSignificance:
                resultList.append([sigFrac-stepSize, mutCnt, oldSig])
                #print('sigFrac={:.2f}, mutCnt={}, sig={:.2f}'.format(sigFrac-stepSize, mutCnt, oldSig))
                break
            oldSig = sig
               
    resultDf = pd.DataFrame(resultList, columns=['sigFrac', 'mutCnt', 'significance'])
    if plotResult:
        resultDf.plot.scatter(x='mutCnt', y='sigFrac')
    return resultDf
    

def findLogRatioFromMinCnt(fermDf, sigVal, minCnt):
    """
    given a min count, a significance and a dataframe with fermentation data, 
    find the desired log ratio that will allow minCnt number of mutations to be significant
    """
    df = sigVsMutCnt(sigVal)
    row = df[df['mutCnt'] == minCnt]['sigFrac']
    if len(row) > 0:
        sigFrac = row.iloc[0]
    
    logRatio = findLogRatio(fermDf, sigFrac)
    return logRatio
    
    
def findLogRatio(fermDf, desiredSigFrac):
    """
    determine the minimum log ratio needed to achieve a certain significance fraction
    fermDf is the output of getAgFermData from get_ferm_data_v4.py
    sigFrac is a float with the desired significance fraction value
    
    For most analyses, we are interested in being able to detect significant
    results with 2, 3 or 4 mutations
    """
    #desiredSigFrac = 0.22
    stepSize = 0.001
    
    minRatio = 0
    maxRatio = 10
    #minRatio = abs(fermDf['logRatio']).min()
    #print(minRatio)
    #maxRatio = abs(fermDf['logRatio']).max()
    #print(maxRatio)
    totalNum = len(fermDf)
    
    for ratio in np.arange(minRatio, maxRatio, stepSize):
        numWithPhenotype = len(fermDf[abs(fermDf['logRatio']) >= ratio])
        sigFrac =  numWithPhenotype/totalNum
        
        if sigFrac < desiredSigFrac:
            print('fraction with phenotype: {}/{}   sigFrac={:.2f}, logRatio={:.2f}'.format(
                    numWithPhenotype, 
                    totalNum, 
                    sigFrac, 
                    ratio))
            break
    return ratio