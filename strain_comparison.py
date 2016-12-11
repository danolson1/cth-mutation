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
    
    # add new columns
    inMut['isOrigin'] = originMutSer # flag to indicate if a mutation is an origin mutation
    
    # the combination of a unique parent strain ID and origin mutation ID indicates the appearance
    # of a mutation.  Sometimes several strains have both the same origin mutation ID and same parent ID
    # this usually means that the mutation occurred once, and that the strains are sisters
    inMut['mut-par'] = inMut['mutID'].astype(str) + '-' + inMut['ParentID'] 
    
    return inMut


def cleanOriginMut(originMutDf, inMut, strainTbl):
    """
    start with origin mutant dataframe
    find repeated origin mutants and delete them
    return the resulting originMut dataframe
    """
    badMut = findRepeatOriginMut(originMutDf, inMut, strainTbl)
    badMutList = badMut['mutID'].tolist()
    
    result = originMutDf.loc[~originMutDf['mutID'].isin(badMutList), :]
    print('cleanOriginMut: repeated origin mutations removed')
    
    return result


def findRepeatOriginMut(originMutDf, inMut, strainTbl):
    """
    find origin mutations that are also present earlier in the lineage
    these are likely to be wrong because:
      1. either the mutation is called incorrectly 
      (mutations on the threshold of calling appear and disappear) or
      2. the strain lineage is wrong
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
        while (parId != rootStrainID) and (parId is not np.nan):
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


def getStrainList(inDf):
    """
    get list of unique strains from a dataframe of mutations
    the list of strains is either in a column with the name
    'Strain' or 'StrainID' 
    """
    strCol = findStrainColumn(inDf)
    strainList = inDf[strCol].unique().tolist()
    return strainList
 
    
def findStrainColumn(inDf):
    """ find the column in a dataframe that has the strain ID numbers """
    if 'Strain' in inDf.columns:
        result = 'Strain'
    elif 'StrainID' in inDf.columns:
        result = 'StrainID'
    else:
        raise ValueError ('No Strain or StrainID columns found')
    
    return result
    

def getStrainLineage(strainList, strainTbl):
    """
    calculate the lineages for all of the strains in the strain list
    """
    strainCol = findStrainColumn(strainTbl)
    
    parChiDf = findParentList(strainList, strainTbl)
    
    lineageSer = pd.Series([False] * len(parChiDf)) 
    
    for index, row in parChiDf.iterrows():
        parList = []
        parId = row['ParentID']
        #print(row['StrainID'])
        while parId is not None: 
            parList.append(parId)
            # check to see if a parent strain was found
            parResult = parChiDf.loc[parChiDf[strainCol] == parId, 'ParentID']
            if len(parResult) > 0: # if the series is not empty, take the first value
                parId = parResult.iloc[0]
            else:
                parId = None
    
        lineageSer.loc[index] = parList
    
    parChiDf['Lineage'] = lineageSer 

    return parChiDf


def getOneStrainLineage(strainID, strainTbl, inMut = None):
    """ 
    get the lineage for a single strain 
    strainId is the strain whose lineage is requested
    inMut contains a list of strains to be considered for lineage construction (if None, use all strains)
    strainTbl (i.e. LLStrain.xlsx) has the master list of parent-child relationships
    """
    
    # build strain list
    if inMut is not None:
        strainList = getStrainList(inMut)
        #print('inMut used')
    else:
        strainList = getStrainList(strainTbl)
        #print('strainTbl used')
    
    sl = getStrainLineage(strainList, strainTbl)
    result = sl.loc[sl['StrainID'] == strainID, 'Lineage'].iloc[0]
    
    return result
    

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


def make_distance_matrix(allAnno, minReadFrac = 0.80):
    """
    starting with a dataframe of annotated mutations ('Chromosome', 'readFrac', 'mutID', 'Strain'),
    make 3 matricies
    distMat is a distance matrix between each strain (Hamming distance)
    newMutMat is a matrix with the number of new mutations in each 'child' strain compared to its 'parent'
    lostMutMat is a matrix with the number of mutations each 'child' strain has lost compared to each 'parent'
    minReadFrac is set to filter out marginal reads.  Default value is 0.80
    """
    # list of strains with zero distance (i.e. identical strains)
    zeroDistList = []
    
    # make list of strains, need to append LL1004 because it doesn't have any mutations, by definition 
    strList = allAnno['Strain'].unique()
    # if LL1004 isn't in the list, add it
    if 'LL1004' not in strList:
        strList = np.append(strList, 'LL1004')
    
    # clean up list of mutations
    rightGenome = allAnno['Chromosome'] == 'Cth_DSM_1313_genome'
    rightReadFrac = allAnno['readFrac'] >= minReadFrac
    # use this for finding 'new' mutations, false positives are less of a problem, 
    # and some adaptation experiments give mutations with low penetration (i.e. readFrac < 0.9)
    cl = allAnno.loc[rightGenome, :] 
    
    # use this for finding 'lost' mutations, 
    # include additional filter by readFraction to avoid false positives
    cl2 = cl.loc[rightReadFrac, :] 
    
    # make empty distance matrix
    # these are global variables
    #global distMat
    #global newMutMat
    #global lostMutMat
    distMat = pd.DataFrame(None, index=strList, columns = strList )
    newMutMat = pd.DataFrame(None, index=strList, columns = strList )
    lostMutMat = pd.DataFrame(None, index=strList, columns = strList )
    newMutID = pd.DataFrame(None, index=strList, columns = strList )
    lostMutID = pd.DataFrame(None, index=strList, columns = strList )
    
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
            
            newMutSet = chiSetAll.difference(parSetAll)
            newMutID.set_value(par, chi, newMutSet)
            newMutMat.set_value(par, chi, len(newMutSet))
            
            lostMutSet = parSetHigh.difference(chiSetHigh)
            lostMutID.set_value(par, chi, lostMutSet)
            lostMutMat.set_value(par, chi, len(lostMutSet))
            
            distVal = len(newMutSet) + len(parSetAll.difference(chiSetAll) )
            distMat.set_value(par, chi, distVal )
            
    return (distMat, newMutMat, lostMutMat, newMutID, lostMutID)
            
    
def findZeroDistStrainPairs(distMat):
    """
    find pairs of strains that share all mutations in common
    these are likely to be 'sister' colonies
    """
    df = pd.DataFrame(distMat.index, columns = ['Strain'])
    df['joinCol'] = 1 # temporary column for merging
    df2 = pd.merge(df, df, on='joinCol')
    df2['dist'] = df2.apply(lambda row: distMat.loc[row['Strain_x'], row['Strain_y']], axis=1)
    
    sameStrain = df2['Strain_x'] == df2['Strain_y']
    goodDist = df2['dist'] == 0
    df3 = df2.loc[~sameStrain & goodDist, :]
    df3 = df3.copy()
    
    # make sets of left-right pairs to allow us to get rid of duplicates (regardless of the order)
    df3['set'] = df3.apply(lambda row: frozenset([row['Strain_x'], row['Strain_y']]), axis=1)
    df3.reset_index(drop = True, inplace = True)
    df4 = df3.loc[:,'set']
    df4.drop_duplicates(inplace = True)
    df5 = df3.loc[df3.index.isin(df4.index), ['Strain_x', 'Strain_y']]
    result  = df5
    
    return result
    
    
def compare_mutations(parStr, chiStr, mutDf, mutType):
    """
    given a parent strain (parStr)
    a child strain (chiStr)
    and a dataframe with mutations (mutDf)
    and the type of comparison (mutType) mutType can be either 'lost' or 'new'
    make a list of the 'lost' or 'new' mutations for manual analysis"""
    
    parMutSet = frozenset(mutDf.loc[mutDf['Strain'] == parStr, 'mutID'].tolist())
    chiMutSet = frozenset(mutDf.loc[mutDf['Strain'] == chiStr, 'mutID'].tolist())
    
    # list of columns to include at the end
    resultColList = ['mutID', 'Strain', 'Chromosome', 'Source', 
                    'StartReg', 'EndReg', 'Description', 'Type', 
                     'Annotation name', 'readFrac']
    
    if mutType == 'lost':
        mutList = parMutSet.difference(chiMutSet)
        inMutList = mutDf['mutID'].isin(mutList)
        rightStr = mutDf['Strain'] == parStr
    elif mutType == 'new':
        mutList = chiMutSet.difference(parMutSet)
        inMutList = mutDf['mutID'].isin(mutList)
        rightStr = mutDf['Strain'] == chiStr
    else:
        print('error, mutType not specified correctly')

    # mutations present in the parent strain that were lost in the child strain 
    print('# differences = ', len(mutList))

    result = mutDf.loc[inMutList & rightStr, resultColList].sort_values('StartReg')
    result['readFrac'] = result['readFrac'].round(2) # clean up decimal places
    return result    
 
    
def findBestParentStrain(chiStrain, distMat, newMutMat, lostMutMat, newID, lostID):
    """
    given a strain, find the strain most likely to be its parent, based on mutation data
    algorithm summary:
      1. find the set of strains with the fewest number of lost mutations
      2. in that set, find the strain with the fewest number of new mutations
      3. note that ties are possible
    return a dataframe with all of the possible matches sorted from best to worst
    """
    resultList = [] # empty list to hold candidate parent strains
    chiColumn = lostMutMat.loc[lostMutMat.index != chiStrain, chiStrain] # prevent strain from searching for itself as a parent
    minLost = chiColumn.sort_values().unique() # list of numbers of lost mutations
    # find strains that match the minLost criteria
    counter = -1 # start at -1 because the first pass through the loop always adds 1 
    #print('minLost= ', minLost)
    for i in minLost:
        #print('i= ', i)
        minNew = -1 # minimum number of new mutations, -1 indicates counter has been reset
        strIdx = chiColumn[chiColumn == i].index
        #print('strIdx ', strIdx)
        # make a submatrix from the newMutations matrix
        # choose the strain with the fewest new mutations
        m = newMutMat.loc[strIdx, chiStrain].sort_values()
        #print('-----------')
        #display(m)
        for index, row in m.iteritems():
            #print(index, counter)
            if row > minNew:
                #resultList.append([index, counter, i, row, distMat.loc[index, chiStrain]])
                minNew = row
                counter += 1
            
            resultList.append([index, counter, i, row, distMat.loc[index, chiStrain], lostID.loc[index, chiStrain], newID.loc[index, chiStrain]])
    
    result = pd.DataFrame(resultList, columns=['Strain', 'Order', 'num lost', 'num new', 'Distance', 'lostMutID', 'newMutID']).set_index('Strain')
    
    return result
    
    
def findBestParentStrainV2(chiStrain, distMat, newMutMat, lostMutMat, newID, lostID):
    """
    given a strain, find the strain most likely to be its parent, based on mutation data
    algorithm summary:
      1. find the set of strains with shortest Hamming distance
      2. go down the list until you find a strain where numNew > numLost
      3. that is the best parent strain
      
      2. in that set, find the strain with the fewest number of new mutations
      3. note that ties are possible
    return a dataframe with all of the possible matches sorted from best to worst
    """
    resultList = [] # empty list to hold candidate parent strains
    chiColumn = lostMutMat.loc[lostMutMat.index != chiStrain, chiStrain] # prevent strain from searching for itself as a parent
    minLost = chiColumn.sort_values().unique() # list of numbers of lost mutations
    # find strains that match the minLost criteria
    counter = -1 # start at -1 because the first pass through the loop always adds 1 
    #print('minLost= ', minLost)
    for i in minLost:
        #print('i= ', i)
        minNew = -1 # minimum number of new mutations, -1 indicates counter has been reset
        strIdx = chiColumn[chiColumn == i].index
        #print('strIdx ', strIdx)
        # make a submatrix from the newMutations matrix
        # choose the strain with the fewest new mutations
        m = newMutMat.loc[strIdx, chiStrain].sort_values()
        #print('-----------')
        #display(m)
        for index, row in m.iteritems():
            #print(index, counter)
            if row > minNew:
                #resultList.append([index, counter, i, row, distMat.loc[index, chiStrain]])
                minNew = row
                counter += 1
            
            resultList.append([index, counter, i, row, distMat.loc[index, chiStrain], lostID.loc[index, chiStrain], newID.loc[index, chiStrain]])
    
    result = pd.DataFrame(resultList, columns=['Strain', 'Order', 'num lost', 'num new', 'Distance', 'lostMutID', 'newMutID']).set_index('Strain')
    
    return result
    
        
def buildParentStrainMatrix(distMat, newMutMat, lostMutMat):
    """
    find the best parent strain for each strain
    """
    strainList = distMat.index.tolist()    
    resultMatrix = pd.DataFrame(np.nan, index=distMat.index, columns=distMat.index )
    for strain in strainList:
        parDf = findBestParentStrain(strain, distMat, newMutMat, lostMutMat)
        parSer = parDf['Order'].copy()
        resultMatrix.loc[:, strain] = parSer
        
    return resultMatrix
    
def pairsFromParStrainMat(m):
    """
    given a parent strain matrix from buildParentStrainMatrix
    find the parent strains with score = 0 (i.e. the best ones)
    note, sometimes a strain has more than one parent
    """
    strainList = m.index.tolist()
    resultList = []
    for strain in strainList:
        parSer = m.loc[:,strain]
        bestPar = parSer[parSer == 0].index.tolist()
        resultList.append([strain, bestPar])
    
    result = pd.DataFrame(resultList, columns = ('Strain', 'Best parent'))
    return result
    