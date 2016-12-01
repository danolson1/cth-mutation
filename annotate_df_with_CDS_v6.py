""" 
#####################################################################################################
This module is for determining the coding sequences (CDS) that could be affected by a given mutation
coding sequences can affect 
Upstream mutations could affect gene expression
Mutations within a gene could affect protein function (do we care about silent mutations?)
A given mutation could affect several genes 
a multi-gene deletion
a mutation upstream of 2 genes on opposite strands
a mutation in a location with overalpping genes (on opposite strands)

Inputs:
     - genbank file with CDS annotations (the same file used by CLC for the mutation detection, to
       ensure that the coordinates are correct)
     - maximum upstream distance (integer), probably either 500 bp or 1000 bp
     - dataframe with start and end coordinates of a mutation ('StartReg' and 'EndReg')
     
Outputs:
     - dataframe that lists the mutations and each CDS that might be affected (possibly more than
       one per mutation) 

Dan Olson
10-24-2016
#####################################################################################################
"""

import pandas as pd
from Bio import SearchIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# initialize constants

#uniqueMutDf = pd.read_csv('unique_cth_mutations_10242016.csv')
max_upstream_dist = 500 # maximum upstream distance for a mutation to affect a gene

# get NCBI accession number for chromosomes based on CLC name
# for match_alleleTbl function
chromosomeDict = {'Cth_DSM_1313_genome':'NC_017304.1'}

#--------------------------------------------------------------------------------------
# annotate dataframe based on CDS
#--------------------------------------------------------------------------------------
def annotate_from_gb(inMutDf, genBankFileName, annotateBest = True):
    """
    given a dataframe with mutations and a genbank filename, determine which 
    genes are affected by a given mutation.  If annotateBest is True, for each
    mutation, only choose the best annotation.  If annotateBest is False, for each
    mutation, choose all CDSs that might interact
    """
    
    # first, make the annotation table
    cdsDf = make_annotation_table(genBankFileName)
    
    # then, make a list of unique mutations (to avoid annotating things multiple times)
    uniqueMutDf = find_unique_mutations(inMutDf)
    
    # then use it to annotate the list of mutations
    annoAll = annotate_all_CDS(uniqueMutDf, cdsDf)

    # choose the one best annotation for each mutation
    if annotateBest: 
        result = annotate_best_CDS(uniqueMutDf, annoAll)
    # OR annotate each mutation with all of the potential genes that could be affected
    else: 
        result = annoAll
    
    # map annotations back to original mutations    
    result = map_anno_to_mut(result, inMutDf)
    
    # keep only the specified columns
    if annotateBest:
        # note, the Locus description column can't be added
        # because there may be more than one locus per mutation
        result = result[['mutID',
                         'Strain',
                         'Sample',
                         'Chromosome',
                         'Source',
                         'StartReg',
                         'EndReg',
                         'Description',
                         'Type',
                        #'lBpId',
                        #'rBpId',
                        'Annotation name',
                        'readFrac'
                        ]]
    else:
        result = result[['mutID',
                         'Strain',
                         'Sample',
                         'Chromosome',
                         'Source',
                         'StartReg',
                         'EndReg',
                         'Description',
                         'Type',
                         'Label',
                         'inCDS',
                         'Locus description',
                         'upstreamDist'
                         'readFrac'                        
                        #'lBpId',
                        #'rBpId',
                        #'Locus start',
                        ]]
    
    return result


def match_alleleTbl(uniqueMutDf, alleleTblFilename = 'Allele.xlsx'):
    """ 
    match mutations (unique or all) against known mutations from the allele table 
    return the new dataframe with matched mutations
    also return a dataframe of unmatched mutations that are likely to be 
    targeted modifications
    """
    alleleTbl = pd.read_excel(alleleTblFilename)
    
    # rename the 'Chromosome' in uniqueMutTable to match standard NCBI names
    # currently the 'Chromosome' name is the nickname I use in CLC
    uniqueMutDf['Chromosome'] = uniqueMutDf['Chromosome'].replace(chromosomeDict)
    
    # get just the columns we care about
    alleleSlice = alleleTbl.loc[:, ['ChromosomeID', 'Start', 'End', 'Type', 'Description','AlleleID', 'AlleleName']]
    
    # match mutations
    result = uniqueMutDf.merge(alleleSlice, how='left', left_on = ['Chromosome', 'StartReg', 'EndReg', 'Type', 'Description'], 
               right_on = ['ChromosomeID', 'Start', 'End', 'Type', 'Description'])
    
    # make a table of alleles that weren't matched
    nullAlleleID = result['AlleleID'].isnull()
    notSnp = result['Source'] != 'snp_data'
    notTransp = result['Source'] != 'match_insertion_Cth_t'
    notTdup = result['Type'] != 'Tandem duplication'
    notMatched = result.loc[nullAlleleID & notSnp & notTransp & notTdup, :]
    
    return (result, notMatched)

     

#--------------------------------------------------------------------------------------
# make dataframe with CDS annotations
# **********note, need to make sure chromosome matches, otherwise we'll get wrong annotations
#--------------------------------------------------------------------------------------
def make_annotation_table(genbankFileName):
    """ given a file name, make a table of annotations """
    record  = SeqIO.read(genbankFileName, 'genbank')
    rowList = [] # list to hold rows
    for feature in record.features:
        if feature.type == 'CDS':
            #check to see if the feature has a label
            if 'label' in feature.qualifiers.keys(): 
                l = feature.qualifiers['label'][0]
            else:
                l = ''   
            #check to see if the feature has a locus tag
            if 'locus_tag' in feature.qualifiers.keys(): 
                lt = feature.qualifiers['locus_tag'][0]
            else:
                lt = ''
            
            #check to see if the feature has a note
            if 'note' in feature.qualifiers.keys(): 
                noteStr = feature.qualifiers['note'][0]
                # split into individual strings based on semicolon
                noteList = [x.strip() for x in noteStr.split(';')]
                # split again by colon
                noteList2 = [x.split(': ') for x in noteList]
                noteDict = {x[0]: x[-1] for x in noteList2}
                #print('****noteDict****')
                #for key, value in noteDict.items():
                #    print(key, '\t\t -->', value)
                # choose one key as description
                # first PFAM, if it exists, then TIGRFAM, then KEGG
                if 'PFAM' in noteDict.keys():
                    desc = noteDict['PFAM']
                elif 'TIGRFAM' in noteDict.keys():
                    desc = noteDict['TIGRFAM']
                elif 'KEGG' in noteDict.keys():
                    desc = noteDict['KEGG']
                elif 'SMART' in noteDict.keys():
                    desc = noteDict['SMART']
                else:
                    # if the note parsing didn't work, include the whole string
                    desc = noteStr 
                    print('ANNOTATION PROBLEM the following text in locus %s could not be parsed:\n\"%s\"' %(lt, desc))
            
            
            #store results in a dictionary
            d = {'Label': l,
                 'Locus Tag': lt,
                 'Strand': feature.location.strand,
                 'Start': feature.location.start.position,
                 'End': feature.location.end.position,
                 'Description': desc}
            rowList.append(d)
    cdsDf = pd.DataFrame(rowList)
    # organize columns
    cdsDf = cdsDf[['Label', 'Locus Tag', 'Strand', 'Start', 'End', 'Description']]
    
    # add columns
    cdsDf['Chromosome'] = record.name
    return cdsDf


def find_unique_mutations(inMut):
    """ given a list of mutations from several strains, find the ones that are unique """
    # get a list of the unique mutations from the list of all mutations
    # 11-22-2016 note, previously I included 'Source' (i.e. which part of my python script
    # identified the mutation) for identifying unique mutations, but 
    # I don't actually think that's an important distinction, since the other parts are what 
    # define a unique mutation
    result = inMut.loc[:, ('Chromosome', 'StartReg', 'EndReg', 'Type', 'Description', 'Source')]
    result.drop_duplicates(inplace = True)
    
    # check to see if there are multiple rows that differ only by source
    dupRows = result[result.duplicated(['Chromosome', 'StartReg', 'EndReg', 'Type', 'Description'])]
    if len(dupRows) > 0:
        print("find_unique_mutations ERROR: same mutation from multiple sources identified")
        print(dupRows)
    
    result.sort_values('StartReg', inplace=True)
    result.reset_index(drop = True, inplace = True) # get rid of original index
    result.reset_index(inplace = True) # move new index to a column called 'index'
    result.rename(columns={'index':'mutID'}, inplace=True)
    return result
    

#--------------------------------------------------------------------------------------
# find overlapping mutations
#--------------------------------------------------------------------------------------
def annotate_all_CDS(inUniqueMutDf, inCdsDf):
    """cdsDf
    given a list of annotations and list of unique mutations, find the overlap.
    unique mutations generated by find_unique_mutations function
    for a given mutation, there may be more than one annotation that overlaps
    mutations that are not in or near a CDS are eliminated 
    """
    uniqueMutDf = inUniqueMutDf.copy()
    cdsDf = inCdsDf.copy()
    
    ## create a table with all pairs of mutation and CDS
    # first, add temporary columns to force outer join to create all combinations of rows from each table
    uniqueMutDf['joinIndex'] = 1
    #uniqueMutDf['mutID'] = uniqueMutDf.index.values #already present in CSV file
    cdsDf['joinIndex'] = 1
    mutCdsJoinDf = pd.merge(uniqueMutDf, cdsDf, on = ['joinIndex', 'Chromosome'] , how='inner')
    mutCdsJoinDf.drop('joinIndex', axis=1, inplace=True) # get rid of temporary columns
    
    # determine overlapping genes
    mutCdsJoinDf['fwdUpstreamDist'] = mutCdsJoinDf['Start'] - mutCdsJoinDf['EndReg']
    mutCdsJoinDf['revUpstreamDist'] = mutCdsJoinDf['StartReg'] - mutCdsJoinDf['End']
    mutCdsJoinDf['startInCds'] = ((mutCdsJoinDf['StartReg'] >= mutCdsJoinDf['Start']) 
                                  & (mutCdsJoinDf['StartReg'] <= mutCdsJoinDf['End']))
    mutCdsJoinDf['endInCds'] = ((mutCdsJoinDf['EndReg'] >= mutCdsJoinDf['Start']) 
                                & (mutCdsJoinDf['EndReg'] <= mutCdsJoinDf['End']))
    mutCdsJoinDf['coversCds'] = ((mutCdsJoinDf['StartReg'] <= mutCdsJoinDf['Start']) 
                                 & (mutCdsJoinDf['EndReg'] >= mutCdsJoinDf['End']))
                                 
    # inCds is the list where a mutation is within or covers a CDS
    withinCds = mutCdsJoinDf.loc[(mutCdsJoinDf['coversCds'] | mutCdsJoinDf['startInCds'] | mutCdsJoinDf['endInCds']), :]
    withinCds = withinCds.drop(['fwdUpstreamDist', 'revUpstreamDist'], axis=1)
    
    # select mutations that could be affecting a gene on the forward strand
    upFwdDf = mutCdsJoinDf.loc[(mutCdsJoinDf['Strand'] == 1) 
                                  & (mutCdsJoinDf['fwdUpstreamDist'] >= 0) 
                                  & (mutCdsJoinDf['fwdUpstreamDist'] <= max_upstream_dist), :]
    upFwdDf = upFwdDf.rename(columns={'fwdUpstreamDist':'upstreamDist'})
    upFwdDf.drop('revUpstreamDist', axis=1, inplace=True)
    
    # select mutations that could be affecting a gene on the reverse strand
    upRevDf = mutCdsJoinDf.loc[(mutCdsJoinDf['Strand'] == -1) 
                                  & (mutCdsJoinDf['revUpstreamDist'] >= 0) 
                                  & (mutCdsJoinDf['revUpstreamDist'] <= max_upstream_dist), :]
    upRevDf = upRevDf.rename(columns={'revUpstreamDist':'upstreamDist'})
    upRevDf.drop('fwdUpstreamDist', axis=1, inplace=True)
    
    # add columns to identify sub-tables
    upFwdDf['Source'] = 'upFwdDf'
    upRevDf['Source'] = 'upRevDf'
    withinCds['Source'] = 'withinCds'
    # re-combine tables
    cdsOverlap = pd.concat([upFwdDf, upRevDf, withinCds]).sort_values('Start')
    cdsOverlap.rename(columns={'Description_x':'Description'}, inplace=True)
    
    # select desired columns
    result = cdsOverlap.loc[:, [ 'Chromosome', 
                                 'mutID',           # unique ID for each mutation
                                 'StartReg',
                                 'EndReg',
                                 'Type',            # SNP, SV, etc.
                                 'Description',     # description of the mutation
                                 'Label',           # Gene name (if any), adhE, for example
                                 'Locus Tag',       # Clo1313 number
                                 'Start',           # Starting coordinate of the locus tag
                                 'Source',          # how the mutation is related to the locus (within or upstream of)
                                 'Description_y',   # description of the locus
                                 'upstreamDist'     # distance upstream of locus (if it is not in a locus)
                                ]]
                     
    # rename columns as needed
    result.rename(columns = {'Source':'inCDS',
                             'Description_y': 'Locus description',
                             'Start': 'Locus start'}, inplace=True)
    return result
    

def map_anno_to_mut(annoUniqueMut, inMutDf):
    """ 
    given a set of annotated mutations, map those annotations back
    to the initial dataframe with mutations
    annoUniqueMut can come from annotate_all_CDS
    OR it can come from annotate_best_CDS
    """
    ######### WORKING ON THIS ###############
    result = inMutDf.merge(annoUniqueMut, how='left', on=['Chromosome','StartReg', 'EndReg', 'Type', 'Description'], suffixes = ('_l','_r'))
    
     # rename columns as needed
    result.rename(columns = {'Source_l':'Mut source',
                             'Source_r':'inCDS',
                             'Description_y': 'Locus description',
                             'Start': 'Locus start'}, inplace=True)
    return result


## need to list the columns that are expected
def annotate_best_CDS(uniqueMutDf, allAnnoDf):
    """ 
    given a list of annotated mutations, find the best annotation for each one 
     - uniqueMutDf is a dataframe of unique mutations generated by functions in
        the file process_clc_files_v5.py
     - allAnnoDf a dataframe with all of the CDSs that could be associated with
        any given mutation.  It is generated by annotate_all_CDS
    """
    # make series to hold "annotation name," default value is 'no CDS match'
    annNameSer = pd.Series(data=['no CDS match']*len(uniqueMutDf), index=uniqueMutDf.index )
    
    # run annotate_all_CDS to get list of CDS annotations
    #allAnnoDf = adc5.annotate_from_gb(inMutDf, 'Cth DSM 1313 genome-111816v2.gbk')
    #uniqueMutDf = find_unique_mutations(inMutDf)
    
    # loop through unique mutations
    for index, mut in uniqueMutDf.iterrows():
        # make a dataframe for all of the annotations that match the current mutation
        # sometimes nothing matches because the mutation is not near or in a CDS
        annoCds = allAnnoDf.loc[allAnnoDf['mutID'] == index, ('Chromosome', 'Label', 
                                                              'Locus Tag', 'inCDS', 
                                                              'upstreamDist', 'Locus start',
                                                              'mutID')]
        
        # split annoCds into the annotations that are within a CDS and those that are upstream
        annoInCds = annoCds.loc[annoCds['inCDS'] == 'withinCds', :]
        annoUpCds = annoCds.loc[(annoCds['inCDS'] == 'upRevDf') | (annoCds['inCDS'] == 'upFwdDf'), :] 
        
        # if there's only one annotation, and it's in a CDS, that's the best annotation for that mutation
        if len(annoInCds) == 1:
            annNameSer[index] = bestName(annoInCds.iloc[0])
        
        # if the mutation covers multiple CDSs, note the first and last CDS
        elif len(annoInCds) > 1:
            annoInCds.sort_values('Locus start')
            startName = bestName(annoInCds.iloc[0]) # first row
            endName = bestName(annoInCds.iloc[-1]) # last row
            annNameSer[index] = startName + '-' + endName
    
        # if the mutation is not in a CDS, len(annoInCds) == 0
        # see if it's upstream of a CDS
        elif len(annoInCds) == 0 and len(annoUpCds) > 0:
            minDist = int(annoUpCds['upstreamDist'].min())
            
            # check to see if there are multiple CDSs with the exact same upstream distance
            df = annoUpCds.loc[annoUpCds['upstreamDist'] == minDist, :]
            if len(df) > 1: # this case hasn't been tested because it didn't occur in my data set D.O. 11/20/2016
                startName = bestName(df.iloc[0]) # first row
                endName = bestName(df.iloc[-1]) # last row
                annNameSer[index] = str(minDist) + ' bp upstream of ' + startName + ' OR ' + endName
            else:
                annNameSer[index] = str(minDist) + ' bp upstream of ' + bestName(df.iloc[0])
                #raise StopIteration
    
    # add 'annotation name' series onto original mutation dataframe
    uniqueMutDf['Annotation name'] = annNameSer
       
    # choose columns we want to keep
    result = uniqueMutDf[['Chromosome', 
                          'StartReg', 
                          'EndReg', 
                          'Type', 
                          'Description',
                          'mutID', 
                          'Annotation name'
                          ]]
    
    return result
    
def bestName(inRow):
    """ 
    helper function for annotate_best_CDS
    given a row with 'Label' and 'Locus tag', choose the best annotation
    choose 'Label' if it's present
    otherwise choose 'Locus tag'
    return a string with the best name
    """
    if inRow['Label'] == '':
        result = inRow['Locus Tag']
    else:
        result = inRow['Label']
        
    return result
    
    
def output_pivot_table(annoMutDf, strainList, outputFileName):
    """ for a given set of strains, make a pivot table to compare mutations """
       
    # figure out which columns of mutDf to use
    indexCols = ['mutID', 'Chromosome','Source', 'StartReg', 'Description', 'Type', 'Annotation name']
    valueCol = 'readFrac'
    columnCol = 'Strain'
    # combine columns into one list for selectoin
    allColList = indexCols
    allColList.extend(valueCol)
    allColList.extend(columnCol)
    
    # select rows corresponding to desired strains
    mutDf = annoMutDf.loc[annoMutDf['Strain'].isin(strainList), allColList]
    
    # make a pivot table
    p = pd.pivot_table(allMutAnno, 
               index = indexCols, 
               values = valueCol,
               columns = columnCol,
               aggfunc='mean',
               fill_value = 0).round(2) # rounds 'readFrac' to 2 decimal places
               
    # the resulting pivot table is a multi-level index that's hard to work with, resetting the index fixes this
    p.reset_index(inplace = True, drop=True)
    #p.to_excel(outputFileName)
    return p
