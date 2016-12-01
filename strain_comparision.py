""" 
#####################################################################################################
This module is for comparing mutations between strains and determining strain lineages

* find identical strains
* find new mutations (mutations that weren't in the parent strain)
* find 'lost' mutations (mutations that were in the parent strain, but are lost in the child strain)


determining the coding sequences (CDS) that could be affected by a given mutation
coding sequences can affect 
Upstream mutations could affect gene expression
Mutations within a gene could affect protein function (do we care about silent mutations?)
A given mutation could affect several genes 
a multi-gene deletion
a mutation upstream of 2 genes on opposite strands
a mutation in a location with overalpping genes (on opposite strands)

Inputs:
     - dataframe with annotated mutations (output from annotate_df_with_CDS)

     
Outputs:
     - various

Dan Olson
11-30-2016
#####################################################################################################
"""

import os
import pandas as pd
import numpy as np
import sys
sys.path.insert(0, r'C:\Users\Dan\Documents\GitHub\cth-mutation') # for loading other modules from my local github directory
#import annotate_df_with_CDS_v6 as adc6

## constants
# master strain list path
mslPath = r'C:\Users\Dan\Documents\Lynd Lab research\Ctherm CBP project\Big database project  - SFRE\--MANUSCRIPT--\Master strain list 11-29-2016.xlsx'

# distance matricies
# module-level variables because we only expect to have one set per instance
distMat = None
newMutMat = None
lostMutMat = None

def make_distance_matrix(annoDf):
    """
    starting with a dataframe of annotated mutations ('Chromosome', 'readFrac', 'mutID', 'Strain'),
    make 3 matricies
    distMat is a distance matrix between each strain (Hamming distance)
    newMutMat is a matrix with the number of new mutations in each 'child' strain compared to its 'parent'
    lostMutMat is a matrix with the number of mutations each 'child' strain has lost compared to each 'paremt'
    """
    # list of strains with zero distance (i.e. identical strains)
    zeroDistList = []
    
    # make list of strains, need to append LL1004 because it doesn't have any mutations, by definition 
    strList = annoDf['Strain'].unique()
    strList = np.append(strList, 'LL1004') # might be a good idea to check first to see if this is present
    
    # clean up list of mutations
    rightGenome = annoAll['Chromosome'] == 'Cth_DSM_1313_genome'
    rightReadFrac = annoAll['readFrac'] > 0.90
    # use this for finding 'new' mutations, false positives are less of a problem, 
    # and some adaptation experiments give mutations with low penetration (i.e. readFrac < 0.9)
    cl = annoAll.loc[rightGenome, :] 
    
    # use this for finding 'lost' mutations, to avoid false positives
    cl2 = cl.loc[rightReadFrac, :] 
    
    # make empty distance matrix
    # these are global variables
    global distMat = pd.DataFrame(None, index=strList, columns = strList )
    global newMutMat = pd.DataFrame(None, index=strList, columns = strList )
    global lostMutMat = pd.DataFrame(None, index=strList, columns = strList )
    
    # make distance matrix
    for par in strList:
        #print('PAR= ', par)
        for chi in strList:
            # make sets of mutations from cl list (i.e.all)
            parSetAll = frozenset(cl.loc[cl['Strain'] == par, 'mutID'].tolist())
            chiSetAll = frozenset(cl.loc[cl['Strain'] == chi, 'mutID'].tolist())
            
            # make sets of mutations from cl2 list (i.e. only high readFrac values)
            parSetHigh = frozenset(cl2.loc[cl2['Strain'] == par, 'mutID'].tolist())
            chiSetHigh = frozenset(cl2.loc[cl2['Strain'] == chi, 'mutID'].tolist())
            
            
            newMutMat.set_value(par, chi, len(chiSetAll.difference(parSetAll)))
            lostMutMat.set_value(par, chi, len(parSetHigh.difference(chiSetHigh)))
            distVal = len(chiSetAll.difference(parSetAll)) +len(parSetAll.difference(chiSetAll) )
            distMat.set_value(par, chi, distVal )
            
            if (distVal == 0) and (par != chi):
                zeroDistList.append([par, chi].sort()) # sort [par, chi] list because order doesn't matter and it makes removing duplicates easier
                #print('***** zero distance strain pair found *****')
                #print('parent= ', par, ' child= ', chi)
                
    zeroDistDf = pd.DataFrame(zeroDistList, columns = ['Strain1', 'Strain2']).drop_duplicates()