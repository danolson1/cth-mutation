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

#----------------------------------------------------------------------------------------------------
# constants and globals
#----------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------
# functions
#----------------------------------------------------------------------------------------------------

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