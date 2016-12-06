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

# for finding origin mutations
rootStrainID = 'LL1004'

# distance matricies
# module-level variables because we only expect to have one set per instance
distMat = None
newMutMat = None
lostMutMat = None

    
def find_origin_mutations(inMut, strainTbl):
    """
    strainTbl is a table with all of the parent-child relationships, 
    i.e. complete LLStrains table
    it needs to have the following columns:
        "StrainID" with the strain ID and
        "Parent strain" with the parent strain ID
    inMut is a list of mutations.  It needs to have the following columns:
        "mutID," unique mutation identifier
        "Strain," strain ID
    """
    # first find the list of strains that have mutations in the mutation table 
    strainList = inMut['Strain'].unique().tolist()
    
    # get the nearest parent with mutation data
    # the resulting dataframe has 2 columns: StrainID and ParentID
    parChiDf = findParentList(strainList, strainTbl)
    
    # add parent strain column
    inMut = inMut.merge(parChiDf, how = 'left', left_on = 'Strain', right_on = 'StrainID')
    inMut.drop('StrainID', axis=1, inplace=True) # don't need to have this column shown twice
    
    # make series to hold origin mutation flag, start with all false values
    originMutSer = pd.Series([False] * len(inMut)) 
    for index, row in inMut.iterrows():
        # if mutation is not in parent strain, it's an origin mutation
        #print(index, 'checking mutID=', row['mutID'])
        pID = row['ParentID']
        if row['mutID'] not in inMut.loc[inMut['Strain'] == pID, 'mutID'].get_values():
            # found origin mutation
            #print('\tfound origin mutation, mutID=', row['mutID'], ' strain=', row['StrainID'])
            originMutSer.loc[index] = True
    
    inMut['isOrigin'] = originMutSer
    
    return inMut

def findRepeatOriginMut(originMutDf, inMut, strainTbl):
    """
    find origin mutations that are also present earlier in the lineage
    these are likely to be wrong
    either the mutation is called incorrectly (mutations on the threshold of calling appear and disappear)
    or the strain lineage is wrong
    originMutDf is a dataframe with just origin mutations (i.e. output of find_origin_mutations where isOrigin is true)
    inMut is a dataframe with all mutations
    strainTbl (i.e. LLStrains)
    """
    strainList = inMut['Strain'].unique().tolist()
    parChiDf = findParentList(strainList, strainTbl) # parent child list for just the strains that have origin mutations
    
    repeatOriginMuts = []
    
    # loop through origin mutations
    for index, row in originMutDf.iterrows():
        parId = row['ParentID']
        while parId != rootStrainID:
            parentOriginMutList = originMutDf.loc[originMutDf['Strain'] == parId, 'mutID']
            #print(row['mutID'])
            #print(row['Strain'])
            #print(parId)
    
            if row['mutID'] in parentOriginMutList.tolist():
                repeatOriginMuts.append(index)
                print('findRepeatOriginMut ERROR: mutation {0} from strain {1} found in parent strain {2}'.format(row['mutID'], row['Strain'], parId))
    
            parId = parChiDf.loc[parChiDf['StrainID'] == parId, 'ParentID'].iloc[0] # find parent of parent
    
    result = originMutDf.loc[originMutDf.index.isin(repeatOriginMuts), :]
    return result


def getStrainList(inMut):
    """ get list of unique strains from a dataframe of mutations """
    strainList = inMut['Strain'].unique().tolist()
    
    return strainList
    


def getStrainLineage(strainList, strainTbl):
    """
    calculate the lineages for all of the strains in the strain list
    """
    parChiDf = findParentList(strainList, strainTbl)

    lineageSer = pd.Series([False] * len(parChiDf)) 
    
    for index, row in parChiDf.iterrows():
        parList = []
        parId = row['ParentID']
        while parId != rootStrainID:
            parList.append(parId)
            parId = parChiDf.loc[parChiDf['StrainID'] == parId, 'ParentID'].iloc[0]
        parList.append(rootStrainID)
        
        lineageSer.loc[index] = parList
        
    parChiDf['Lineage'] = lineageSer 

    return parChiDf

#--------------------------------------------------------------------------------------
# find all parent strains from list
#--------------------------------------------------------------------------------------
def findParentList(strainList, strainTbl):
    """ 
    given a list of strain ID numbers, and a master strain table (i.e. LLStrains),
    get the closest parent strain
    assume that strain IDs are strings
    empty string '' means parent not found 
    """
    
    # check to see if the strainlist has the root strain, and if not, add it
    if rootStrainID in strainList:
        pass
    else:
        strainList.append(rootStrainID)
        
    parentStrList = []
    
    # loop through all of the strains
    for strain in strainList:
        #print('looking for parent of strain ', strain)
        parId = findParent(strain, strainTbl) # find parent strain
        #print('\tparent is: ', parId)
        
        while parId not in strainList: # is parent in list?
            # if not, then find the parent of the parent
            #print(parId, ' not in list')
            parId = findParent(parId, strainTbl) 
            #print('\tfindParentList: parent not found, new parId=', parId)
            # if the parent strain doesnt not exist, then you can't go any further
            if parId == '':
                break
        parentStrList.append([strain, parId])
    result = pd.DataFrame(parentStrList, columns=['StrainID', 'ParentID'])
    # get rid of rows where parent strain is empty
    result.drop(result[result['ParentID'] == ''].index, inplace = True)
    
    return result

#--------------------------------------------------------------------------------------
# find a single parent strain
#--------------------------------------------------------------------------------------
def findParent(id, strainTbl):
    """ 
    Helper function for findParentsList
    Given a strain ID and table of parent-child relationships
    find the parent strain
        "id" is string with strain ID
        "strainTbl" is dataframe with strain ID and parent strain ID (i.e. LLStrains)
            "StrainID" is strain ID
            "ParentStrain" is parent strain ID
    """
    #print('findParent: id in=', id)
    parId = '' # default "not found" value
    
    parResult = strainTbl.loc[strainTbl['StrainID'] == id, 'ParentStrain'] # might be null
    if len(parResult) > 0:
        parId = parResult.iloc[0]
    else:
        pass
        #print('findParent ERROR: ', id, ' no parent found')
       
    # check if the strain is the wild type strain
    if id == rootStrainID:
        parId = 'wt'
        
    return parId







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
    global distMat
    global newMutMat
    global lostMutMat
    distMat = pd.DataFrame(None, index=strList, columns = strList )
    newMutMat = pd.DataFrame(None, index=strList, columns = strList )
    lostMutMat = pd.DataFrame(None, index=strList, columns = strList )
    
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