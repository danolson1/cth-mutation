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
version 7
3-30-2020
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

# from match_alleleTbl, determines cutoff for matching start and end regions
# Note 12-7-2016, this has been replaced by using isSameRegion for matching and should be removed at some point
maxMismatchDistance = 6 

# for isSameRegion function
minOverlapThresh = 0.98 # mutations must be 98% overlapped to be considered identical

# get NCBI accession number for chromosomes based on CLC name
# for match_alleleTbl function
chromosomeDict = {'Cth_DSM_1313_genome':'NC_017304.1'}

#--------------------------------------------------------------------------------------
# annotate dataframe based on CDS
#--------------------------------------------------------------------------------------
def annotate_from_gb(inMutDf, genBankFileName, annotateBest = True, alleleTblFilename = None):
    """
    given a dataframe with mutations and a genbank filename, determine which 
    genes are affected by a given mutation.  
    If annotateBest is True, for each mutation, only choose the best annotation.  
    If annotateBest is False, for each mutation, choose all CDSs that might interact
    When an allele table filename has been provided, try to match unique mutations to the Allele.xlsx data file
    
    inMutDf: a dataframe with mutation information. Note that the 'Chromosome' name has to match the 'Chromosome' name from the genbank file
        or annotations will not be matched together
    genBankFileName: the source of data about annotations, in Genbank format
    annotateBest:
        True - assigns one locus to each mutation
        False - if multiple loci could be affected by the mutation, all of them are included
    alleleTblFilename: Excel file with alleles to match to data. This comes from the Allele table in my MS Access strain database
    """
    
    # first, make the annotation table
    cdsDf = make_annotation_table(genBankFileName)
    
    # then, make a list of unique mutations (to avoid annotating things multiple times)
    uniqueMutDf = find_unique_mutations(inMutDf)
    # and map mutation IDs back to the list of all mutations
    inMutDf = map_unique_mutid_to_mut(uniqueMutDf, inMutDf)
    
    # then annotate the unique list of mutations with CDS overlaps
    # each unique mutID may have several annotations
    uniqueAllAnno = annotate_all_CDS(uniqueMutDf, cdsDf)

    # choose the one best annotation for each mutation
    if annotateBest: 
        uniqueAnno = annotate_best_CDS(uniqueMutDf, uniqueAllAnno) 
        """# clean up columns
        # note, the Locus description column can't be added
        # because there may be more than one locus per mutation
        uniqueAnno = uniqueAnno[['mutID', # allMut
                                 'Strain', # allMut
                                 'Sample', # allMut
                                 'Chromosome', # allMut
                                 'Source', # allMut
                                 'StartReg', # allMut
                                 'EndReg', # allMut
                                 'Description', # allMut
                                 'Type', # allMut
                                #'lBpId',# allMut
                                #'rBpId', # allMut
                                 'Annotation name',
                                 'readFrac' # allMut
                                ]]
        """
    # OR annotate each mutation with all of the potential genes that could be affected
    else: 
        uniqueAnno = uniqueAllAnno
        """
        uniqueAnno = uniqueAnno[['mutID', # allMut
                                 'Strain',# allMut
                                 'Sample',# allMut
                                 'Chromosome',# allMut
                                 'Source',# allMut
                                 'StartReg',# allMut
                                 'EndReg',# allMut
                                 'Description',# allMut
                                 'Type',# allMut
                                 'Label',
                                 'inCDS',
                                 'Locus description',
                                 'upstreamDist'
                                 'readFrac'                        
                                #'lBpId',
                                #'rBpId',
                                #'Locus start',
                               ]]
        """
    
    # map annotations back to original mutations based on the mutation ID    
    inMutAnno = inMutDf.merge(uniqueAnno, how = 'left', on='mutID', suffixes = ('', '_y')) # use _y flag to keep track of duplicate columns
    inMutAnno = inMutAnno.loc[:, ~inMutAnno.columns.str.contains('_y')] # get rid of columns with _y flag
    
    # match Alleles, if the Allele.xlsx location has been specified
    if alleleTblFilename is not None:
        uniqueMutAllele = match_alleleTbl(uniqueMutDf, alleleTblFilename)
        # match AlleleID and AlleleName back to the original mutants
        inMutAnno = inMutAnno.merge(uniqueMutAllele, how = 'left', on='mutID', suffixes = ('', '_y')) # use _y flag to keep track of duplicate columns
        inMutAnno = inMutAnno.loc[:, ~inMutAnno.columns.str.contains('_y')] # get rid of columns with _y flag
    
    result = inMutAnno

    return result


def match_alleleTbl(uniqueMutDf, alleleTblFilename = 'Allele.xlsx'):
    """ 
    match mutations (unique or all) against known mutations from the allele table 
    return only the important columns, mutID, AlleleID and AlleleName
    """
    alleleTbl = pd.read_excel(alleleTblFilename)
    
    # rename the 'Chromosome' in uniqueMutTable to match standard NCBI names
    # currently the 'Chromosome' name is the nickname I use in CLC
    uniqueMutDf['Chromosome'] = uniqueMutDf['Chromosome'].replace(chromosomeDict)
    
    # get just the columns we care about
    alleleSlice = alleleTbl.reindex(['ChromosomeID', 'Start', 'End', 'Type', 'Description','AlleleID', 'AlleleName'], axis = 'columns')
    
    # match mutations
    result = uniqueMutDf.merge(alleleSlice, how='left', left_on = ['Chromosome', 'Type', 'Description'], 
               right_on = ['ChromosomeID', 'Type', 'Description'])
               
    # use isSameRegion to find close matches for start and end regions 
    result['isSameRegion'] = result.apply(lambda x: isSameRegion(x['StartReg'], x['EndReg'], 
                                                                 x['Start'], x['End']), axis=1)
    result = result.loc[result['isSameRegion'], :]
    # drop extra columns added to calculate same region
    result.drop(['Start', 'End', 'isSameRegion'], axis=1, inplace = True)    
    
    
    
    """
        # get rid of other columns from uniqueMut, except the ones we care about
    uniqueMut = uniqueMut.loc[:,['mutID', 'Chromosome', 'Description', 'Type', 'StartReg', 'EndReg'] ]
    
    # if inMutDf already has a column called 'mutID' we need to get rid of it
    if 'mutID' in inMutDf.columns:
        inMutDf.drop('mutID', axis=1, inplace=True)
    
    # merge based on the columns that should match exactly
    result = inMutDf.merge(uniqueMut, how='left', on=['Chromosome', 'Type', 'Description'])
    

    
    # rename columns as needed
    result.rename(columns = {'StartReg_x' : 'StartReg',
                             'EndReg_x' : 'EndReg'}, 
                  inplace = True)
    
    # move the column to head of list using index, pop and insert
    cols = list(result)
    cols.insert(0, cols.pop(cols.index('mutID')))
    result = result.loc[:, cols]
               
    result['sDist'] = result['StartReg'] - result['Start']
    result['eDist'] = result['EndReg'] - result['End']
    result['totDist'] = result['sDist'].abs() + result['eDist'].abs()
    
    closeMatch = result['totDist'] < maxMismatchDistance # allow a total difference of 6 bp for matches
    result = result.loc[closeMatch, :]
    """
    
    # make a table of alleles that weren't matched
    # this should go it its own function
    """
    umdf = uniqueMutDf.copy()
    umdf = umdf.merge(result.loc[:, ['AlleleID', 'mutID']], on='mutID', how='left')
    nullAlleleID = umdf['AlleleID'].isnull()
    
    notSnp = umdf['Source'] != 'snp_data'
    notTransp = umdf['Source'] != 'match_insertion_Cth_t'
    notTdup = umdf['Type'] != 'Tandem duplication'
    rightChrom = umdf['Chromosome'] == chromosomeDict['Cth_DSM_1313_genome']
    notMatched = umdf.loc[nullAlleleID & notSnp & notTransp & notTdup & rightChrom, :]
    
    uniqueNotMatched = notMatched.loc[:, ['mutID', 'Chromosome', 
                                          'StartReg', 'EndReg', 'Description', 
                                          'Type', 'Annotation name']]
    uniqueNotMatched.drop_duplicates(inplace = True)
    """
        
    return result


def map_unique_mutid_to_mut(uniqueMut, inMutDf):
    """ 
    given a set of unique mutation IDs (uniqueMut),
    map them on to a set of mutations (inMut)
    based on exact match of Chromosome, Type and Description
    and close match to StartReg and EndReg
    
    uniqueMut columns:
      "Chromosome"
      "Description"
      "Type"
      "StartReg"
      "EndReg"
      "mutID"
      
    inMutDf columns:
      "Chromosome"
      "Description"
      "Type"
      "StartReg"
      "EndReg"
      "mutID"
      and maybe other columns too
    """
    
    # get rid of other columns from uniqueMut, except the ones we care about
    uniqueMut = uniqueMut.loc[:,['mutID', 'Chromosome', 'Description', 'Type', 'StartReg', 'EndReg'] ]
    
    # if inMutDf already has a column called 'mutID' we need to get rid of it
    if 'mutID' in inMutDf.columns:
        inMutDf.drop('mutID', axis=1, inplace=True)
    
    # merge based on the columns that should match exactly
    result = inMutDf.merge(uniqueMut, how='left', on=['Chromosome', 'Type', 'Description'])
    
    # use isSameRegion to find close matches for start and end regions 
    result['isSameRegion'] = result.apply(lambda x: isSameRegion(x['StartReg_x'], x['EndReg_x'], 
                                                                 x['StartReg_y'], x['EndReg_y']), axis=1)
    result = result.loc[result['isSameRegion'], :]
    # drop extra columns added to calculate same region
    result.drop(['StartReg_y', 'EndReg_y', 'isSameRegion'], axis=1, inplace = True)
    
    # rename columns as needed
    result.rename(columns = {'StartReg_x' : 'StartReg',
                             'EndReg_x' : 'EndReg'}, 
                  inplace = True)
    
    # move the column to head of list using index, pop and insert
    cols = list(result)
    cols.insert(0, cols.pop(cols.index('mutID')))
    result = result.loc[:, cols]

    return result
    

#--------------------------------------------------------------------------------------
# make dataframe with CDS annotations
# **********note, need to make sure chromosome matches, otherwise we'll get wrong annotations
#--------------------------------------------------------------------------------------
def make_annotation_table(genbankFileName):
    """ given a file name, make a table of annotations """
    rowList = [] # list to hold rows
    records = list(SeqIO.parse(genbankFileName, 'genbank')) # to deal with multiple DNA sequences in a single file
    for record in records:
        for feature in record.features:
            if feature.type == 'CDS':
                desc = None # default value  Without this you get an error when you try to store disc in a dictionary later
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
                
                # if the feature has a label, but no description, set the description as the label
                if desc is None:
                    desc = l
                
                #store results in a dictionary
                d = {'Label': l,
                     'Locus Tag': lt,
                     'Strand': feature.location.strand,
                     'Start': feature.location.start.position,
                     'End': feature.location.end.position,
                     'Description': desc,
                     'Chromosome': record.name
                     }
                rowList.append(d)
    cdsDf = pd.DataFrame(rowList)
    # organize columns
    cdsDf = cdsDf[['Chromosome', 'Label', 'Locus Tag', 'Strand', 'Start', 'End', 'Description']]
    
    # add columns
    #cdsDf['Chromosome'] = record.name
    return cdsDf


def isSameRegion(start_x, end_x, start_y,  end_y):
    """
    determine whether 2 mutations are actually identical,
    even if they have slightly different start and end coordinates
    usually this is due to differences in read mapping
    sometimes it can be caused by the presence of two adjacent mutations
    i.e. a SNP and a breakpoint
    """
    isSame = False # boolean flag to hold the output, whether the two mutations are the same
    
    totDiff = abs(start_x - start_y) + abs(end_x - end_y) # total difference between start and end regions
    if totDiff == 0: # if the region is identical, no further analysis needed
        isSame = True
        return isSame 
    else:  # otherwise, check to see how big the difference is 
        len_x = (end_x - start_x)
        len_y = (end_y - start_y)
        maxLen = max(len_x, len_y)
        if maxLen > 0:
            minOverlap = 1 - (float(totDiff) / maxLen)
            if minOverlap > minOverlapThresh:
                isSame = True
    
    return isSame

def find_unique_mutations(inMut):
    """ given a list of mutations from several strains, find the ones that are unique """
    # get a list of the unique mutations from the list of all mutations
    # 11-22-2016 note, previously I included 'Source' (i.e. which part of my python script
    # identified the mutation) for identifying unique mutations, but 
    # I don't actually think that's an important distinction, since the other parts are what 
    # define a unique mutation
    result = inMut.loc[:, ('Chromosome', 'StartReg', 'EndReg', 'Type', 'Description')]
    result.drop_duplicates(inplace = True)
    
    # clean up result
    result.sort_values('StartReg', inplace=True)
    result.reset_index(drop = True, inplace = True) # get rid of original index
    result.reset_index(inplace = True) # move new index to a column called 'index'
    result.rename(columns={'index':'mutID'}, inplace=True)
    
    # check for duplicate mutID, where the StartReg and EndReg regions are close but not identical
    # make list of duplicated mutations
    checkDup = result.merge(result, on=['Chromosome', 'Type', 'Description'])
    checkDup['isSameRegion'] = checkDup.apply(lambda x: isSameRegion(x['StartReg_x'], x['EndReg_x'], x['StartReg_y'], x['EndReg_y']), axis=1)
    notIdentical = checkDup['mutID_x'] != checkDup['mutID_y']
    sameRegion = checkDup['isSameRegion']
    checkDup = checkDup.loc[notIdentical & sameRegion, :]
    dupMutIdx = checkDup['mutID_x'].unique().tolist() # list of duplicated mutation indices
    
    # eliminate those rows
    dupIndList = [] # list of duplicate indices to remove
    df = result.loc[result.index.isin(dupMutIdx), :] # make a dataframe with just the set of potential duplicate rows
    df.set_index('mutID')
    
    for i in range(len(dupMutIdx)):
        test = dupMutIdx[i]
        #print('checking mutation ', test)
        for j in range(len(dupMutIdx)-(i+1)):
            candidate = dupMutIdx[j+i+1]
            #print('\t against mutation ', candidate)
            CTDmatch = df.loc[test, ['Chromosome', 'Type', 'Description']].equals(df.loc[candidate, ['Chromosome', 'Type', 'Description']])
            regionMatch = isSameRegion(df.loc[test,'StartReg'], df.loc[test,'EndReg'], df.loc[candidate,'StartReg'], df.loc[candidate,'EndReg'] )
            
            if (CTDmatch and regionMatch):
                #print('\t\t mutation matched')
                dupIndList.append(candidate)
       
    #print(set(dupIndList))
    result = result.loc[~result['mutID'].isin(dupIndList), :] # list of mutations with no duplicates
    
    
    
    return result

#--------------------------------------------------------------------------------------
# find overlapping mutations
#--------------------------------------------------------------------------------------
def annotate_all_CDS(inUniqueMutDf, inCdsDf):
    """
    given a list of annotations and list of unique mutations, find the overlap.
    unique mutations generated by find_unique_mutations function
    for a given mutation, there may be more than one annotation that overlaps
    mutations that are not in or near a CDS are eliminated 
    
    inCdsDf is a dataframe that lists all of the possible annotations and has the columns:
        'Start'
        'End'
        'Chromosome'
        'Description' (description of the locus)
        
    inUniqueMutDf is a dataframe that lists the mutations that were found and has columns:
        'Chromosome'
        'StartReg'
        'EndReg'
        'Description' (description of the mutation)
    """
    print(inCdsDf.columns)
    print(inUniqueMutDf.columns)
    
    
    uniqueMutDf = inUniqueMutDf.copy()
    cdsDf = inCdsDf.copy()
    
    ## create a table with all pairs of mutation and CDS
    # first, add temporary columns to force outer join to create all combinations of rows from each table
    uniqueMutDf['joinIndex'] = 1
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
    cdsOverlap = pd.concat([upFwdDf, upRevDf, withinCds], sort = True).sort_values('Start')
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

    


### NEED TO DELETE THIS FUNCTION, NO LONGER USED ###
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
    
    # make series to hold "Locus description"
    locusDescriptionSer = pd.Series(data=['']*len(uniqueMutDf), index=uniqueMutDf.index )

    # make series to hold "Best locus tag"
    # this is useful for correlating other data (like gene expression data) with the mutation
    bestLocusTagSer = pd.Series(data=['']*len(uniqueMutDf), index=uniqueMutDf.index )
    
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
                                                              'mutID', 'Locus description')]
        
        # split annoCds into the annotations that are within a CDS and those that are upstream
        annoInCds = annoCds.loc[annoCds['inCDS'] == 'withinCds', :]
        annoUpCds = annoCds.loc[(annoCds['inCDS'] == 'upRevDf') | (annoCds['inCDS'] == 'upFwdDf'), :] 
        
        # if there's only one annotation, and it's in a CDS, that's the best annotation for that mutation
        if len(annoInCds) == 1:
            annNameSer[index] = bestName(annoInCds.iloc[0])
            locusDescriptionSer[index] = annoInCds.iloc[0].loc['Locus description']
            bestLocusTagSer[index] = annoInCds.iloc[0].loc['Locus Tag']
        
        # if the mutation covers multiple CDSs, note the first and last CDS
        elif len(annoInCds) > 1:
            annoInCds.sort_values('Locus start')
            startName = bestName(annoInCds.iloc[0]) # first row
            endName = bestName(annoInCds.iloc[-1]) # last row
            annNameSer[index] = startName + '-' + endName
            locusDescriptionSer[index] = 'multiple loci'
            bestLocusTagSer[index] = annoInCds.iloc[0].loc['Locus Tag'] # use the locus tag of the first row
            
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
                locusDescriptionSer[index] = df.iloc[0].loc['Locus description']
                bestLocusTagSer[index] = df.iloc[0].loc['Locus Tag'] 
                #raise StopIteration
    
    # add 'annotation name' series onto original mutation dataframe
    uniqueMutDf['Annotation name'] = annNameSer
    
    uniqueMutDf['Locus description'] = locusDescriptionSer
    
    uniqueMutDf['Locus Tag'] = bestLocusTagSer # keep this named 'Locus Tag' to maintain compatibility with annotateAll
       
    # choose columns we want to keep
    result = uniqueMutDf[['Chromosome', 
                          'StartReg', 
                          'EndReg', 
                          'Type', 
                          'Description',
                          'mutID', 
                          'Annotation name',
                          'Locus Tag',
                          'Locus description'
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
    
    
def output_pivot_table(annoMutDf, strainList, outputFileName = 'outputPivotTable.xlsx'):
    '''
    # Note 06-23-2017 : this function is currently not working
    # The most up-to-date version of the code is in a Jupyter notebook
    # C:\\Users\\Dan\\Documents\\Lynd Lab research\\Ctherm CBP project\\Hfs operon in T sacch - Aysenur\\
    # resequencing data analysis\Python analysis\hfsB strains reseq analysis 6-22-2017.ipynb
    # End note

    #Note: 'Annotation name' is a column that is produced by the annotateBest function

    for a given set of strains, make a pivot table to compare mutations
    '''
       
    # figure out which columns of mutDf to use
    indexCols = ['mutID', 'Chromosome','Source', 'StartReg', 'Description', 'Type', 'Annotation name']
    valueCol = 'readFrac'
    columnCol = 'Strain'
    # combine columns into one list for selectoin
    allColList = indexCols.copy()
    allColList.append(valueCol)
    allColList.append(columnCol)
    
    # select rows corresponding to desired strains
    mutDf = annoMutDf.loc[annoMutDf['Strain'].isin(strainList), allColList]
    
    # make a pivot table
    p = pd.pivot_table(annoMutDf, 
               index = indexCols, 
               values = valueCol,
               columns = columnCol,
               aggfunc='mean',
               fill_value = 0).round(2) # rounds 'readFrac' to 2 decimal places
               
    # the resulting pivot table is a multi-level index that's hard to work with, resetting the index fixes this
    p.reset_index(inplace = True, drop=True)
    #p.to_excel(outputFileName)
    return p
