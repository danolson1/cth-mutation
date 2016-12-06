""" 
#####################################################################################################
Functions for processing CLC files
    ## get input from a variety of places ##
        * Manesh
        * batch processing of CLC files
    
    ## identify breakpoints by BLAST ##
    
    ## match breakpoints together to identify mutations ##
        * match breakpoints together to identify insertions, deletions and transposons
        * compare breakpoints with SV data (CLC is good at calling some SV mutations)
        * compare deletions with read coverage data (a true deletion should have zero coverage in the deleted region)
        * combine with SNP data
        * final output is csv file with all mutations


Inputs:
     - location of CLC files with mutation information (SNP, BP, coverage, SV, indel)
     - location of fasta file with insertions

Outputs:
     - csv file with all mutations
     - csv file of matched BPs
     - csv file of unmatched BPs

Note on column names:
NAME	    DESCRIPTION

Strain	    LL strain number (or other strain identifier)
Sample	    Name of DNA sequence data file from JGI
StartReg	Starting coordinate of mutation, based on Cth genome
EndReg	    Ending coordinate of mutation, based on Cth genome
lBpId	    Left breakpoint ID number
rBpId	    Right breakpoint ID number
readFrac	Fraction of reads supporting the presence of the mutation (float from 0-1)
Type	    Type of mutation (deletion, insertion, targeted deletion, etc.)
Description	Name of insertion sequence.  This can be empty and may be filled in later
Source	    Python function that generated the result
AlleleID	Used for matching known mutations to expected mutations to confirm strain construction

Dan Olson
11-16-2016
version 4
Uses Anaconda3 as the python interpreter
#####################################################################################################
"""


# Imports
import os
import re
import subprocess # for run_command function
import pandas as pd

import time
from Bio import SearchIO
from Bio import SeqIO
from Bio.Alphabet import IUPAC   # the alphabet we want is IUPACUnambiguousDNA
from Bio.Blast import NCBIWWW
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

#initialize global parameters
##########################################################
# constants

# maximim distance between breakpoint associated with a deletion 
# and where region of <10x coverage starts
# originally I had this set to 10, but in AG1155, I needed to increase it to 12
maxDeletionDistanceError = 20

maxHitEndDistance = 20 # maximum distance between BLAST hit and end of sequence, important for insertions
 
eValue = 1e-6 # default E-value to use for BLAST searches
minNumReads = 4 # for filtering breakpoints
minLength = 15 # for filtering breakpoints

# folder with output text files from CLC SFRE batch run
clcDataFilePath = r'F:\Seq data for automated analysis\CLC output files'

# filenames for output results
identifiedBpFilename = 'identified_breakpoints.csv'
unidentifiedBpFilename = 'unidentified_breakpoints.csv'
matchedBpFilename = 'matched_breakpoints.csv'
unmatchedBpFilename = 'unmatched_breakpoints.csv'

# default blastSubjectList to use if we don't specify a different one
# nr gets appended automatically by make_blast_subject_list
blastSubjectFiles = ['insertion_sequences.fa', # pyk(Tsc), pta merodiploid and Shuen's pathway insertions
                    'Cth_transposons.fa',
                    'Cth_known_regions.fa', 
                    'Cth_homology.fa', 
                    'Cth_1313_CDS_annotations.fa', 
                    'Cth_DSM_1313_genome.fa']

# constants for identifying real mutations from threshold data
# see cleanLowThreshold function
minLowCoverageLength = 50 # ignore low coverage regions shorter than this
maxAverageReadCoverage = 1 # ignore regions with maximum coverage higher than this

############################################################


def process_files(clcFilePath, identifiedBp = None):
    """
    main method for this file
    get mutation data from CLC text files
    then analyze mutations
    return dataframe with combined mutations
    
    Option: supply a dataframe with identified breakpoints to avoid having to re-search
    """
    # read CLC text files, store output in a dictionary of results
    resultDict = get_clc_data_from_txt(clcFilePath)
    
    
    inBpDf = resultDict['BP']
    svDf = resultDict['SV']
    lowThresh = resultDict['Thr_low10']
    highThresh = resultDict['Thr_high20']
    snpDf = resultDict['SNP']
    indelDf = resultDict['InDel']
    
    result = combineMutations( inBpDf = inBpDf,
                               svDf = svDf,
                               lowThresh = lowThresh,
                               highThresh = highThresh,
                               snpDf = snpDf,
                               indelDf = indelDf,
                               identifiedBp = identifiedBp)
    return result
    

def combineMutations(inBpDf = None, svDf = None, lowThresh = None, highThresh = None, snpDf = None, indelDf = None, identifiedBp = None):
    """ 
    idBpDf: dataframe with identified breakpoints (from identify_breakpoints_v2.py)
    svDf: dataframe with structural variant mutations from CLC
    lowThresh: dataframe with areas of less than 10x coverage, from CLC
    snpDf: dataframe with SNP variants from CLC
    indelDf: dataframe with indel variants from CLC
    identifiedBp: breakpoints that have already been identified by BLAST, 
         this is useful if you've already done this step 
         and don't want to re-do it, since it takes a long time
    all inputs are optional 
    """
    # initialize variables
    cleanSv = None
    cleanBp = None
    cleanSnp = None
    
    # process breakpoint data first
    if inBpDf is not None:
       cleanBp = cleanBreakpoints(inBpDf)
       # make BLAST subject list
       bsl = make_blast_subject_list(blastSubjectFiles, addNr = True)
       # search for breakpoints only if they have not been provided 
       if identifiedBp is None:
           (identifiedBp, unidentifiedBp) = nameBreakpointWithBlast(breakPointDf = cleanBp,
                                                                    blastSubjectList = bsl ) 
    ## for troubleshooting when I don't want to re-run the BLAST search
    #identifiedBp = inBpDf # for troubleshooting
    #cleanBp = 1 # for troubleshooting
    
    # clean up input data
    # check to make sure it's present first
    print('\n### cleaning input data ###')
    if (svDf is not None) and (lowThresh is not None):
        cleanSv = cleanStructuralVariants(svDf, lowThresh) # make sure deletions have < 10 coverage across most of their length
    if (snpDf is not None):
        cleanSnp = cleanSnpMutations(snpDf)
    if (indelDf is not None):
        cleanI = cleanIndel(indelDf)
    if (highThresh is not None):
        (cleanHighThresh, avgCov) = cleanHighThreshold(highThresh)

    # match breakpoints
    # previously I had been eliminating breakpoints after I matched them
    # this created some problems where the order of searching was important
    # so I changed the search function to always search with all of the breakpoints
    # potentially this can give a problem where the same breakpoint is assigned to 
    # multiple mutations
    
    # match transposon insertions 
    print('\n### matching %s breakpoints to transposons ###' %(str(len(identifiedBp))))
    (matchedIns, remainingBp2) = match_insertions(identifiedBp, 
                                                  maxDist = 50, 
                                                  seqList = 'Cth_transposons.fa', 
                                                  filterByDist = True, 
                                                  filterByEndDist = False)   
    print('\n### matching %s breakpoints to insertions ###' %(str(len(identifiedBp))))
    # match other insertions, i.e. merodiploid DNA
    (matchedIns2, remainingBp3) = match_insertions(identifiedBp, 
                                                   maxDist = 5000, 
                                                   seqList = 'insertion_sequences.fa', 
                                                   filterByDist = True, 
                                                   filterByEndDist = True)   
 
    # check to make sure cleanSv and cleanBp are both present
    if (cleanSv is not None) and (cleanBp is not None):
        # match breakpoints to clean SV mutations, keep only SVs that have matched breakpoints, 
        # collect unmatched breakpoints
        print('\n### matching %s breakpoints to sv ###' %(str(len(identifiedBp)))) 
        (matchedSv, remainingBp1) = match_bp_to_sv(cleanSv, identifiedBp)
          
        # keep the SV mutations who didn't match any breakpoints
        # the evidence that they're real is less strong than for the matchedSv mutations, 
        # but we may still be interested in them
        leftoverCleanSv = cleanSv.loc[~cleanSv.index.isin(matchedSv.index), ['Strain', 'Sample', 'Chromosome', 'StartReg', 'EndReg', 'Variant ratio', 'Name', 'Source']]
        leftoverCleanSv.rename(columns = {'Name':'Type', 'Variant ratio':'readFrac'}, inplace = True)
        
        # check for duplicates based on close distance
        leftoverCleanSv = checkForDuplicate(matchedSv, leftoverCleanSv)
        
           

    # export list of unmatched breakpoints
    # need to fix this part, since I got rid of the hierarchical searching
    #remainingBp3.to_csv('unmatched_breakpoints.csv')
    
    # combine mutations from all dataframes (except low threshold)
    # need to check for duplicate mutations, strain LL1380 is good test case for this
    result = pd.concat([cleanSnp, cleanI, matchedIns, matchedIns2, matchedSv, leftoverCleanSv, cleanHighThresh])

    # check cleanLowThresh for duplicates
    if (lowThresh is not None):
        cleanLowThresh = cleanLowThreshold(lowThresh)
        cleanLowThresh = checkForDuplicate(result, cleanLowThresh)
        result = pd.concat([result, cleanLowThresh])
     
    # match to known mutations and identify with unique mutation ID
    
    # export final dataframe as "all_cth_mutations"
    result.reset_index(inplace=True, drop=True)
    return result
    

def checkForDuplicate(inGood, inPossibleDupe):
    """ 
    given 2 dataframes, check for possible duplicate mutations
    match chromosome, strain, sample, type
    StartReg and EndReg of both duplications should be within 50 bp
    return a dataframe with the non-duplicated rows
    """
    result = inPossibleDupe.copy().reset_index()
    result = result.merge(inGood, how='inner', on=['Chromosome', 'Sample', 'Strain', 'Type'])
    result['sDist'] = result['StartReg_x'] - result['StartReg_y']
    result['eDist'] = result['EndReg_x'] - result['EndReg_y']
    result['totDist'] = result['sDist'].abs() + result['eDist'].abs()
    result.sort_values('totDist', inplace = True) # sort values so that next step takes the smallest result
    result = result.groupby('index').first() # choose the rows with the shortest totDist
    result[['Chromosome', 'StartReg_x', 'StartReg_y', 'totDist']]
    
    # if the total distance between the two start regions and the two end regions is less than 20
    # the two regions probably have a substantial amount of overlap
    dupRows = result['totDist'] < 20
    
    # get just the good rows
    goodRows = result.loc[~dupRows, :]
    
    # return the non-duplicated rows from the original inPossibleDupe dataframe
    result = inPossibleDupe.loc[inPossibleDupe.index.isin(goodRows.index), :]
    return result
    
    
def cleanIndel(rawIndel):
    """ 
    get rid of heterozygous rows
    eliminate unneeded columns 
    split region into start and end regions
    fix chromosome label
    
    """
    result = rawIndel.loc[rawIndel['Zygosity'] == 'Homozygous', ['Strain',
                                                                'Sample',
                                                                'Chromosome', 
                                                                'Region', 
                                                                'Type', 
                                                                'Variant ratio']]
    
    # split Region into StartReg and EndReg
    result = splitRegion(result)
    result.drop('Region', axis=1, inplace=True)
    
    # fix chromosome spaces
    result['Chromosome'] = result['Chromosome'].apply(fixChromosomeSpaces)
    
    #result.rename(columns = {
     #       }, inplace=True)
    
    # make new Source column
    result['Source'] = 'indel_data'
    
    return result
    
    
def cleanBreakpoints(rawBP):
    """
    given a dataframe with breakpoints, get rid of the ones that don't meet certain criteria 
        filter the breakpoints to speed up the BLAST search
        minLength: don't bother doing a blast search if BP is shorter than this
        minNumReads: don't bother doing a blast search if there are fewer than this number of reads
                   initially I was using 10, but came across a case where an important BP only
                   had 9 reads.
    """
    print('cleaning up identifiedBpDf...')
    
    rawBP = splitRegion(rawBP)
    # fix chromosome spaces
    rawBP['Chromosome'] = rawBP['Chromosome'].apply(fixChromosomeSpaces)
    
    # make sure BP sequence (in 'Unaligned' column) only contains the letters A, G, C and T
    rawBP['Unaligned'] = rawBP['Unaligned'].replace({'[^AGCT]':''}, regex=True)
    
    # build boolean filters
    seqNotNull = pd.notnull(rawBP['Unaligned'])  # check to make sure the 'Unaligned' value is not null
    seqLongEnough = rawBP['Unaligned'].str.len() > minLength  # choose a minimum BP length,
                                                           # short sequences don't have enough information for BLAST
                                                           # seq length should be determined from length of 
                                                           # 'Unaligned', not the 'Unaligned length' integer
                                                           # since this may have changed based on regex replacements
    enoughReads = rawBP['Reads'] > minNumReads  # samples with only a few reads are likely to be sequencing errors
    cleanRawBP = rawBP.loc[seqNotNull & seqLongEnough & enoughReads, :]

    print('done cleaning up rawBP')
    return cleanRawBP
    
        
def splitRegion(inDf):
    """ 
    helper method for splitting a CLC "Region" field into "StartReg" and "EndReg" parts 
    input a dataframe with one column called Region 
    returns a dataframe with new columns: "StartReg" and "EndReg" 
    """
    pattern = re.compile('(\d*)') # find a group of digits of any length
    startSer = pd.Series([None]*len(inDf), index=inDf.index)
    endSer = pd.Series([None]*len(inDf), index=inDf.index)
    for index, row in inDf.iterrows():
        regStr = str(row['Region'])
        #print(regStr)
        # find groups of digits in each 'Region' string
        # note, sometimes 'Region' is an integer
        # this usually returns 2 groups of digits, but sometimes 4 (for join-type regions)
        # the output includes a bunch of empty strings as well
        result = pattern.findall(regStr)
        result = list(filter(None, result)) # eliminate empty strings
        startSer[index] = int(result[0]) # call the first regex match the start 
        #print('\t', result[0])
        # check to make sure the result list isn't empty
        if len(result) >= 2: # if there are 2 ore more regex matches, call the 2nd one the end
            i = 1
        if len(result) <= 1: # if there's only one regex match, use that one for both start and end
            i = 0
        endSer[index] = int(result[i]) # call the second regex match the end
        #print('\t', result[i])
        

    # add start and end regions as columns to rawSv dataframe
    inDf['StartReg'] = startSer
    inDf['EndReg'] = endSer
    return inDf


def cleanHighThreshold(highThresh):
    """ 
    clean up high threshold data.  This can be used to identify the presence of exogenous DNA such as repB and cat
    also calculates the average read depth of the C. therm genome
    """
    highThresh = splitRegion(highThresh)
    highThresh['Distance'] = highThresh['EndReg'] - highThresh['StartReg']
    
    # estimate average coverage of C. therm chromosome
    # for regions that have 20x or higher coverage 
    # (i.e. ignore deletion regions)
    cthChrom = highThresh['Chromosome'] == 'Cth DSM 1313 genome'
    cthThresh = highThresh.loc[cthChrom, :]
    sumOfDistance = cthThresh['Distance'].sum()
    avgCoverage = (cthThresh['Average value in window'] * (cthThresh['Distance'] / sumOfDistance)).sum()
    
    # calculate copy number for all regions
    highThresh['CopyNumber'] = highThresh['Average value in window'] / avgCoverage
    
    notCthChrom = highThresh['Chromosome'] != 'Cth DSM 1313 genome'
    goodDistance = highThresh['Distance'] > minLowCoverageLength
    result = highThresh.loc[notCthChrom & goodDistance, ['Strain', 'Sample', 'Chromosome', 'StartReg', 'EndReg', 'CopyNumber']]

    # fix chromosome spaces
    result['Chromosome'] = result['Chromosome'].apply(fixChromosomeSpaces)

    result['Source'] = 'high_threshold'
    
    return (result, avgCoverage)

    
def cleanLowThreshold(lowThresh):
    """ clean up low threshold data.  This can be used to identify deletions due to lack of sequence coverage """
    # split start and end regions
    lowThresh = splitRegion(lowThresh)
    lowThresh['Distance'] = lowThresh['EndReg'] - lowThresh['StartReg']
    
    # distance needs to be long enough
    goodDistance = lowThresh['Distance'] > minLowCoverageLength # constant value set at top of file
    
    # only search the C. therm chromosome
    rightChrom = lowThresh['Chromosome'] == 'Cth DSM 1313 genome'
    
    # real deletion regions should have an average coverage of less than 0.1
    lowReadCoverage = lowThresh['Average value in window'] < maxAverageReadCoverage # constant value set at top of file
    
    # since average coverage is 50-100x, 1-readCoverage is a reasonable approximation of the readFraction
    lowThresh['readFrac'] = 1-(lowThresh['Average value in window'] / 50)
    
    result = lowThresh.loc[goodDistance & rightChrom & lowReadCoverage, ['Strain', 'Sample', 'Chromosome', 'StartReg', 'EndReg', 'readFrac']]
    result['Source'] = 'low_threshold'
    
    # fix chromosome spaces
    result['Chromosome'] = result['Chromosome'].apply(fixChromosomeSpaces)
    
    result['Type'] = 'Deletion'
    
    return result

    
def cleanStructuralVariants(rawSv, lowThresh):
    """ clean up structural variants
    split Region into Start and End
    keep only mutations that are "deletion" or "tandem duplication" types
    make sure that deletions have zero coverage in the deletion region 

    #### NOTE NEED TO TEST WITH TANDEM DUPLICATION ####
    """
    print('cleaning structural variations...')
    
    # split Region into start and end
    rawSv = splitRegion(rawSv)
    lowThresh = splitRegion(lowThresh)
    
    # get rid of rows from the rawSv data with a split group number
    # split groups indicate complex sv results, and seem to be wrong, in general
    rawSv = rawSv.loc[rawSv['Split group'].isnull(), :]
    
    # get rid of rows, except for the following types
    # where Name = Deletion, Evidence = Cross mapped breakpoints
    # where Name = Insertion, Evidence = Tandem duplication
    # note: need to use &, | instead of and, or,  see http://stackoverflow.com/questions/8632033
    isDeletion = (rawSv['Name'] == 'Deletion') & (rawSv['Evidence'] == 'Cross mapped breakpoints')
    isInsertion = (rawSv['Name'] == 'Insertion') & (rawSv['Evidence'] == 'Tandem duplication')
    isReplacement = (rawSv['Name'] == 'Replacement')
    rawSv = rawSv.loc[isDeletion | isInsertion | isReplacement, :]
    
    # check the read coverage of deletions, for each strain
    strainList = rawSv['Strain'].unique()
    correctDeletionList = []
    badDeletionList = []
    
    # loop through strains
    for strain in strainList:
        #print('************ working on strain %s ... **************************'%(strain))
        tempStrainSv = rawSv.loc[rawSv['Strain'] == strain, :]
        tempStrainCov = lowThresh.loc[lowThresh['Strain'] == strain, :]
        chromosomeList = tempStrainSv['Chromosome'].unique()
        
        # loop through each chromosome within each strain
        for chrom in chromosomeList:
            #print('    working on chromosome %s ...'%(chrom))
            tempSv = tempStrainSv.loc[tempStrainSv['Chromosome'] == chrom, :]
            tempCov = tempStrainCov.loc[tempStrainCov['Chromosome'] == chrom, :]
        
            # limit dataframe to SV mutations that are deletions
            tempSvDel = tempSv.loc[isDeletion, :]
            
            # make a dataframe with the best match for each SV row
            tsd = tempSvDel.loc[:,['Strain', 'StartReg', 'EndReg']].reset_index() # temporary SV dataframe
            tc = tempCov.loc[:,['Strain', 'StartReg', 'EndReg']].reset_index() # temporary coverage dataframe
            t2 = tsd.merge(tc, on='Strain', suffixes=('_SV', '_cov')) # match all SV rows to all coverage rows
            t2['startDiff'] = t2['StartReg_cov'] - t2['StartReg_SV'] # distance between start of SV and start of low coverage region
            t2['endDiff'] = t2['EndReg_SV'] - t2['EndReg_cov'] # distance between end of low coverage region and end of SV
            t2['totDiff'] = t2['startDiff'].abs() + t2['endDiff'].abs() # total for sorting
            bestMatchDf = t2.sort_values('totDiff', ascending = True).groupby('index_SV').first()
            
            # split this dataframe based on whether it meets or does not meet the maxDeletionDistanceError threshold
            # note that we only need a one-sided test
            # it's fine if the low coverage region extends beyond the SV region (i.e. when startDiff or endDiff is negative)
            startIsGood = bestMatchDf['startDiff'] <= maxDeletionDistanceError
            endIsGood = bestMatchDf['endDiff'] <= maxDeletionDistanceError
            
            bestMatchDfGood = bestMatchDf.loc[(startIsGood & endIsGood), :]
            bestMatchDfBad = bestMatchDf.loc[~(startIsGood & endIsGood), :]
            
            # add the index values to the correct list
            correctDeletionList.extend(bestMatchDfGood.index.tolist())
            badDeletionList.extend(bestMatchDfBad.index.tolist())
           
    # make new dataframe with correct deletions
    correctDeletions = rawSv.loc[correctDeletionList, :].copy()
    correctDeletions['Source'] = 'sv_data low_coverage'
    
    badDeletions = rawSv.loc[badDeletionList, :].copy()
    badDeletions['Source'] = 'sv_data'
    
    # make dataframe with insertions
    insertionDf = rawSv.loc[isInsertion, :].copy()
    insertionDf['Source'] = 'sv_data'
    
    # originally I was planning to eliminate the 'badDeletions' because 
    # some of them show up due to transposon inversions
    # however this signature also shows up in strains where adhE was deleted and then re-inserted
    # i.e. LL1153
    result = pd.concat([correctDeletions, badDeletions, insertionDf])
    
    # after cleaning, reset the row numbering to be consecutive
    result.reset_index(drop = True, inplace = True) 
    
    # fix chromosome spaces
    result['Chromosome'] = result['Chromosome'].apply(fixChromosomeSpaces)
    
    return result


def cleanSnpMutations(rawSnp):
    """ clean SNP data """

    # split 'Region' field into start and end
    rawSnp = splitRegion(rawSnp)
    
    rawSnp['readFrac'] = rawSnp['Frequency']*(1/100) # convert 0-100 scale of Freqency to 0-1 scale of readFrac
    rawSnp['Amino acid change'].replace('.*p\.', '', regex = True, inplace = True) # clean up amino acid change text
    
    # make new Source column
    rawSnp['Source'] = 'snp_data'
    
    # combine 'Reference' 'Allele' and 'Amino acid change in the description field
    # all fields have a reference and allele
    rawSnp['Description'] = rawSnp['Reference'] + ' --> ' + rawSnp['Allele'] 
    
    # split the list into rows that have an 'Amino acid change' value
    hasAaChange = ~rawSnp['Amino acid change'].isnull()
    a = rawSnp.loc[hasAaChange, 'Description'] + ', ' + rawSnp.loc[hasAaChange, 'Amino acid change'] # rows with amino
                                                                                                     # acid change
    b = rawSnp.loc[~hasAaChange, 'Description'] # rows without an amino acid change
    c = a.append(b, verify_integrity = True).sort_index()
    rawSnp['Description'] = c
    
    # fix chromosome spaces
    rawSnp['Chromosome'] = rawSnp['Chromosome'].apply(fixChromosomeSpaces)
    
    # select columns to keep
    colList = [  'Strain',
                 'Sample',
                 'Chromosome',
                 'StartReg',
                 'EndReg',
                 'Type',
                 'Description',
                 'readFrac',
                 'Source']
    cleanRawSnp = rawSnp[colList]
    return cleanRawSnp    
    

def match_insertions(inDf, maxDist, seqList, filterByDist = False, filterByEndDist = False):
    # select just the breakpoints with BLAST hits from the list of sequences from seqName
    # for example, if seqName is 'Cth_transposons.fa', just match transposons
    transpDf = inDf[inDf['BLAST source'] == seqList]
    
    # check to see if transpDf is empty
    if len(transpDf) == 0:
        print('No breakpoints supplied to match_insertions. '
              'This will generate a runtime error when you try to merge 2 empty dataframes. '
              'This error can be ignored')
    
    
    # split transposons into left and right groups and only keep the relevant columns
    transpDfLeft = transpDf.loc[transpDf['Name'] == 'Left breakpoint', ['Strain','Sample', 'Chromosome', 
                                                                        'Hit Name', 'Seq Name', 'StartReg', 
                                                                        'EndReg', 'Fraction non-perfectly mapped', 
                                                                        'BP distance', 'Hit end distance']]
    transpDfRight = transpDf.loc[transpDf['Name'] == 'Right breakpoint', ['Strain','Sample','Chromosome', 
                                                                          'Hit Name', 'Seq Name', 'StartReg',
                                                                          'EndReg', 'Fraction non-perfectly mapped', 
                                                                          'BP distance', 'Hit end distance']]
    
    
    # do an inner join of left and right sets of breakpoints
    # match all breakpoints with the same name and chromosome
    lr = pd.merge(transpDfLeft, transpDfRight, how='inner', on=['Strain','Sample','Chromosome', 'Hit Name'], suffixes=(' L', ' R'))
    
    # filter by distance
    # choose left and right breakpoint pairs based on how close they are
    # this is useful for transposons, 
    # its only important for targetd mutations if there are multiple identical insertions in the same strain
    if filterByDist:
    # transposons have a small negative distance
        lr['absDist'] =  lr['EndReg R'] - lr['StartReg L']
        lr['absDist'] = lr['absDist'].abs()
        correctDist = (lr['absDist'] < maxDist)
        lr = lr.loc[correctDist, :]
        # if there are multiple matches, try to choose the best one (or at least get rid of some bad ones)
        # select all rows that match the smallest distance
        minDistIdx = lr.groupby('Seq Name L')['absDist'].transform(min) == lr['absDist']
        lr = lr.loc[minDistIdx, :]    
        
    # filter by end distance (i.e. the distance between where the breakpoint starts matching and 
    # the start of the sequence).  Real insertion sequences have breakpoints that match only at the
    # ends (i.e. Hit end distance is close to 0).  This is not important for other types of mutations
    if filterByEndDist:
        correctRDist = (lr['Hit end distance R'].abs() < maxHitEndDistance)
        correctLDist = (lr['Hit end distance L'].abs() < maxHitEndDistance)
        lr = lr.loc[correctRDist & correctLDist, :]
        # if there are several matches for a given breakpoint, choose the one with the shortest end distances
        # allEndDist is a composite of the left and right BP distances and hit end distances
        # a low score means that the breakpoint starts matching the target sequence close to the chromosomal junction and close to the 
        # end of the insertion sequence
        lr['allEndDist'] = lr['BP distance L'].abs() + lr['BP distance R'].abs() + lr['Hit end distance L'].abs() + lr['Hit end distance R'].abs()
        # select rows where 'allEndDist' is equal to the minimum 'allEndDist' value in the group
        minAllDistIdx = lr.groupby('Seq Name L')['allEndDist'].transform(min) == lr['allEndDist']
        lr = lr.loc[minAllDistIdx, :]

    result = lr.copy()
    ## Clean up dataframe
    # add new columns as necessary
    result['readFrac'] = (result['Fraction non-perfectly mapped L'] + result['Fraction non-perfectly mapped R']) / 2
    result['Type'] = 'Insertion'
    result['Source'] = 'match_insertion_' + seqList[:5] # add first 6 letters of seq list to distinguish origins 
    
    # rename columns
    # for deletions, the StartReg coordinate will always be lower than the EndReg coordinate
    # for insertions, the StartReg coordinate will be greater than the EndReg coordinate
    result.rename(columns = {'Hit Name':'Description',
                             'Seq Name L':'lBpId',
                             'Seq Name R':'rBpId',
                             'EndReg R': 'StartReg',
                             'StartReg L':'EndReg'
                            }, inplace=True)
    
    
    # fix chromosome spaces
    result['Chromosome'] = result['Chromosome'].apply(fixChromosomeSpaces)
    
    # select only the desired columns for the final output
    result = result[['Strain','Sample','Chromosome','StartReg','EndReg','lBpId','rBpId','readFrac','Type','Description','Source']]
    
    # figure out which breakpoints weren't matched
    # make a list of the values in the lBpId and rBpId columns
    matchedBpList =  sorted(result['lBpId'].astype(int).tolist() + result['rBpId'].astype(int).tolist())
    # match those against the original list
    matchedBooleanSeries = inDf['Seq Name'].astype(int).isin(matchedBpList)
    # find the original breakpoints whose breakpoint IDs are not in the matchedBpList
    remainingBp = inDf.loc[~matchedBooleanSeries, :]
    
    return (result, remainingBp)


def match_bp_to_sv(inSv, inBpDf):
    """
    match breakpoints to SV mutations
    note, at some point we might want to also match these against the table of known mutations
    this would give us an "AlleleID" column
    """
    tsv = inSv.copy()
    tsv['Left breakpoints'] = tsv['Left breakpoints'].astype(int).astype(str) # cast this as a string for subsequent matching
    
    # split breakpoints into left and right matches
    tbp = inBpDf.copy()
    tbp['Start'] = tbp['StartReg'].astype(str) # cast this as a string to match the string values in 'Left breakpoints' in the SV dataframe
    tbp.drop(['StartReg', 'EndReg'], axis=1, inplace=True)
    tbp.sort_values(['Strain', 'Chromosome', 'Start'], inplace=True, ascending=False)
    tbp.set_index(['Strain', 'Chromosome', 'Start'], inplace=True)
    tbpLeft = tbp.loc[tbp['Name'] =='Left breakpoint', :]
    tbpLeft = tbpLeft.drop('Name', axis=1)
    tbpRight = tbp.loc[tbp['Name'] == 'Right breakpoint', :]
    tbpRight = tbpRight.drop('Name', axis=1)
    
    # join left BPs to SVs
    r1 = tsv.merge(tbpLeft, how='inner', 
                   left_on=['Strain','Chromosome', 'Left breakpoints'], 
                   right_index=True,
                  )
    # join right BPs to SVs
    result = r1.merge(tbpRight, how='inner', 
                  left_on=['Strain','Chromosome', 'Right breakpoints'], 
                  right_index = True,
                  suffixes=('_L','_R')
                 )
    
    # find best match for each SV
    # Hit end distance is the distance from the start of where the BP target starts matching
    #     0 indicates a good match, larger numbers might suggest an incorrect match
    result['totHitEndDist'] = result['Hit end distance_L'].abs() + result['Hit end distance_R'].abs()
    result = result.sort_values('totHitEndDist', ascending = True)
    result = result.groupby(result.index).head(1)
    
    #print('cleaning up dataframe...')
    ## Clean up dataframe
    # rename columns
    result.rename(columns={'Name': 'Type',
                           'Variant ratio': 'readFrac',
                           'Seq Name_L':'lBpId',
                           'Seq Name_R':'rBpId',
                          }, inplace=True)
    
    # rename "Insertion" in "Name" to "Tandem duplication"
    for index, row in result.iterrows():
        if row['Type'] == 'Insertion' and row['Evidence'] == 'Tandem duplication':
            result.set_value(index, 'Type', 'Tandem duplication')
    # add column Source
    result['Source'] = result['Source'] + ' match_bp'
    
    # fix chromosome spaces
    result['Chromosome'] = result['Chromosome'].apply(fixChromosomeSpaces)
    
    # select only the columns we need
    result = result[['Strain', 'Sample', 'Chromosome', 'StartReg', 'EndReg', 'lBpId', 'rBpId', 'readFrac',
                     'Type', 'Source']]
    
    #print('figure out which BPs were not matched')
    # figure out which breakpoints weren't matched
    # make a list of the values in the lBpId and rBpId columns
    matchedBpList =  sorted(result['lBpId'].astype(int).tolist() + result['rBpId'].astype(int).tolist())
    # match those against the original list
    matchedBooleanSeries = inBpDf['Seq Name'].astype(int).isin(matchedBpList)
    # find the original breakpoints whose breakpoint IDs are not in the matchedBpList
    remainingBp = inBpDf.loc[~matchedBooleanSeries, :]   
    return (result, remainingBp)


def fixChromosomeSpaces(inStr):
    """ replace spaces with underscores in 'Chromosome' field 
    this is important because the files imported from CLC have
    spaces, but when a genome is exported as a genbank file,
    the spaces are replaced with underscores
    
    In order to do the CDS matching correctly, I need to match
    chromosomes.  In order for the right chromosomes to match,
    they need to have the exact same spelling"""
    inStr = inStr.replace(' ', '_')
    return inStr

""" 
#####################################################################################################
Functions for getting input data


Inputs:
     - list of files to import (use CLC .csv output)
     - list of fasta files for matching against 
            (note, order matters) 
            (also note, no spaces in filenames)
     - minLength: minimum length of a breakpoint sequence to analyze (default is 15)
     - minNumReads: minimum number of reads for doing a breakpoint analysis (default is 6)

     
Outputs:
     - dataframe with identified breakpoints
     - csv files with identified and with unidentified breakpoints

Dan Olson
11-11-2016
version 2
Uses Anaconda3 as the python interpreter
#####################################################################################################
"""

def get_clc_data_from_txt(filePath = clcDataFilePath):
    """
    Note: needs access to master strain table for relationship between strainID and sampleID
    Given a starting directory with outputs from a CLC SFRE workflow:
     - parse file names to get strain name and mutation type (SV, BP, SNP, etc.)
     - make a table of the files
     - check for duplicate rows, if found use highest numbered version
     - for each type of mutation, make one dataframe to hold the data from all strains
     - group these dataframes together in a list
     - return the list
     """
    # look at the files in the input file path (inFilePath) and import them into a dataframe
    # Note, the file output format changed when I upgraded from CLC workbench v8 to v9
    
    # regex pattern to get strain number, the type of CLC output, and the version number
    # backref 1 = strain ID
    # backref 3 = type of CLC output (threshold, BP, InDel, SV or SNP)
    #         note, SNPs are listed as 'Variants, CTRL, OA, AAC'
    # backref 6 = version.  If empty, this is version 0.  If the workflow is run
    #         multiple times, there will be numbers here.  The highest number is the 
    #         latest version
    p = re.compile( '((AG|LL)\d{3,4})'                          # backref 1 = strain ID
                    '.*(Variants|InDel|BP|SV|threshold)[^._]*'  # backref 3 = type
                    '(_low10|_high20)?'                         # backref 4 - type of threshold (optional)
                    '(\.(\d{1,2}))?'                            # backref 6 - version (optional)
                    '\.txt'
                    )
        
    typeDict = {
        'threshold':'Thr',
        'SV':'SV',
        'BP':'BP',
        'InDel':'InDel',
        'Variants':'SNP'
               }
    
    fileList = []
    for file in os.listdir(filePath):
        if file.endswith('.txt'):
            result = p.findall(file)
            #print(file)
            #print(result)
            strainID = result[0][0]
            fileType = typeDict[result[0][2]]
            
            # add low10 or high20 identifier to the end of Thr
            # low10 is the default option
            if fileType == 'Thr': # get low10 or high20 identifier
                if result[0][3] != '':
                    fileType  = fileType + result[0][3]
                else:
                    fileType  = fileType + '_low10'
            
            if result[0][5] is '': # the file doesn't have a value in the fileVers location
                fileVers = 0
            else:
                fileVers = int(result[0][5]) # get the value
    
            row = [file, strainID, fileType, fileVers]
            fileList.append(row)
        
    clcFilesDf = pd.DataFrame(fileList, columns = ['Filename', 'StrainID', 'Type', 'Version'])
    
    ### for testing
    #clcFilesDf.to_excel('clcFilesDf.xlsx')
    
    # for each strain and type, make sure we only select the row with the highest version number
    maxVersionIndex = clcFilesDf.groupby(['StrainID', 'Type'])['Version'].transform(max) == clcFilesDf['Version']
    clcFilesDf = clcFilesDf.loc[maxVersionIndex, :]
    
    # read master strain list excel file
    # use this to look up sample IDs
    mTbl = pd.read_excel('mTbl.xlsx')
    
    # sometimes this has an 'Index' colum
    # this happens if I opened the mTbl.xlsx file in excel and turned it into a table
    # then pandas reads the initial column as a column called 'Index' not the index
    if 'Index' in mTbl.columns:
        mTbl.drop(['Index'], axis=1, inplace=True)
    
    # make a dataframe for each of the 5 input file types
    typeList = ['BP', 'InDel', 'SV', 'Thr_low10', 'Thr_high20', 'SNP']
    resultDict = dict(zip(typeList, [None]*len(typeList)))
    # loop through types
    for item in typeList:
        print('###########################', item, '###########################')
        resultList = [] # list to hold dataframes from each individual file
        fileTypeDf = clcFilesDf.loc[clcFilesDf['Type'] == item, :]
        
        # loop through each row of the new dataframe
        for index, row in fileTypeDf.iterrows():
            # read the file into a dataframe
            #print(row)
            clcFilePath = os.path.join(filePath, row['Filename'])
            df = pd.read_table(clcFilePath)
            df['Strain'] = row['StrainID']
            sampleName = mTbl.loc[mTbl['Strain ID'] == row['StrainID'], 'Filename'].values[0]
            #print(sampleName)
            df['Sample'] = sampleName
            # need to look up Sample ID from strain table and add it
            resultList.append(df)
            
        # append the individual CLC files of one type into a single dataframe
        clcTypeDf = pd.concat(resultList).reset_index(drop=True)
        
        resultDict[item] = clcTypeDf
        
    return resultDict
     


""" 
#####################################################################################################
Functions for identifying breakpoints from CLC analysis of JGI resequencing data


Inputs:
     - list of files to import (use CLC .csv output)
     - list of fasta files for matching against 
            (note, order matters) 
            (also note, no spaces in filenames)
     - minLength: minimum length of a breakpoint sequence to analyze (default is 15)
     - minNumReads: minimum number of reads for doing a breakpoint analysis (default is 6)

     
Outputs:
     - dataframe with identified breakpoints
     - csv files with identified and with unidentified breakpoints

Dan Olson
11-11-2016
version 2
Uses Anaconda3 as the python interpreter
#####################################################################################################
"""
# this function can probably be deleted DO 11-17-16
def identify_bps(bpFileNames, sequenceFileNames):
    """read files from CLC, merge them, blast against list of sequences and nr, output the result as csv files"""
    # read files from CLC
    df = combine_CLC_files(bpFileNames)
    
    # identify with BLAST
    (matchedDf, unmatchedDf) = nameBreakpointWithBlast(df, sequenceFileNames)
    
    # output the result as a filename
    matchedDf.to_csv(identifiedBpFilename)
    unmatchedDf.to_csv(unidentifiedBpFilename)
    
    return matchedDf
    
# this function can probably be deleted DO 11-17-16
def combine_CLC_files(bpFileNames):
    """ 
    combine text files from CLC output into dataframes
    one dataframe for breakpoints, one for SVs, one for SNPs
    """
    
    # read in the files from the bpFileNames list
    bpDfList = []
    p = re.compile('(\d{3,4})') # extract the first string of 3 or 4 digits and call that the strain number
    for fileName in bpFileNames:
        df = pd.read_csv(fileName)
        df['Strain'] = int(p.findall(fileName)[0]) # add the strain name as a column
        bpDfList.append(df)  
    rawBP = pd.concat(bpDfList, ignore_index = True) # concatenate them all together
    return cleanBreakpoints(rawBP) # clean up the output and return

# this function can probably be deleted DO 11-17-16
def importManeshFiles():
    """ import data from Manesh text files, convert to dataframe, concatenate and filter """
    # the raw BP data from Manesh came in 2 files, one for my data and one for Adam Guss's data
    rawBPdan = pd.read_table(r"October2016\DanOlson_SFRE_V5_Test3\Strain_Sample_BrkPnts.Test3.txt.Tabbed")
    rawBPadam = pd.read_table(r"October2016\AdamGuss_SFRE_V5_Test3\Strain_Sample_BrkPnts.Test3.txt.Tabbed")
    

    # concatenate the data from Dan and from Adam in a single frame, rawBP
    frames = [rawBPdan, rawBPadam]
    rawBP = pd.concat(frames, ignore_index=True)

    # filter the breakpoints to speed up the BLAST search
    ##seqNotNull = pd.notnull(rawBP['Unaligned'])  # check to make sure the 'Unaligned' value is not null
    ##seqLongEnough = rawBP['Unaligned length'] > minLength  # choose a minimum BP length,
    ## short sequences don't have enough information for BLAST
    ##enoughReads = rawBP['Reads'] > minNumReads  # samples with only a few reads are likely to be sequencing errors
    ##cleanRawBP = rawBP[seqNotNull & seqLongEnough & enoughReads]
    ##return(cleanRawBP)
    return (cleanBreakpoints(rawBP))


def make_blast_subject_list(local_files, addNr = True):
    """ given a list of fasta filenames, make a blast list """
    # helper for name breakpoints
    # list of BLAST subjects
    # first column is filename (for local) or database name (for NCBI)
    # second column is true for local and false for NCBI
    # all of these files need to be placed in the working directory
    # Example below:
    # blastSubjectList = [['Cth_known_regions.fa', True],
    #                     ['Cth_transposons.fa', True],
    #                     ['Cth_homology.fa', True],
    #                     ['Cth_1313_CDS_annotations.fa', True],
    #                     ['Cth_DSM_1313_genome.fa', True],
    #                     ['nr', False]]
    blastSubjectList = []
    for file in local_files:
        blastSubjectList.append([file, True]) 
    
    if addNr == True:
        blastSubjectList.append(['nr', False]) # add the nr database as the final item in the list
    return blastSubjectList


def nameBreakpointWithBlast(breakPointDf, blastSubjectList):
    """
    assign names to breakpoints by sequential BLAST search against local fasta files and nr
        breakPointDf is a dataframe with the breakpoints to be matched
        blastSubjectList holds the subject sequences and iscreated by the function make_blast_subject_list
    Overview:
    1. breakpoints are searched against each local file in order
    2. matching breakpoints are named
    3. remaining breakpoints are searched against the next sequence
    4. breakpoints remaining after all local_files have been searched are searched against nr
    5. matched and unmatched breakpoints are written as csv files and also returned
    """   
    #make blastSubjectList manually to have control of whether or not to include the nr search
    #blastSubjectList = make_blast_subject_list(local_files) # first make "subject list" from local filenames
    #blastSubjectList = local_files

    remainingBPs = breakPointDf.copy()
    setOfFoundBPs = []
    for [item, localSearch] in blastSubjectList:
        print("\n\n\n###### BLAST search " + str(len(remainingBPs)) + " sequences against " + item + '\t #######')
        # only do a BLAST search if there is one or more sequences in the remainingBPs list
        if len(remainingBPs) >= 1:
            knownDf = localBLASTn(remainingBPs, subjectSeqFileName=item, localSearch=localSearch)

            # sort cleanRawBP into 2 dataframes, matched and unmatched, need to add .copy() to get new frame instead of slice
            matchedDf = remainingBPs.loc[remainingBPs.index.isin(knownDf['Seq Name'])].copy()
            unmatchedDf = remainingBPs.loc[~remainingBPs.index.isin(knownDf['Seq Name'])].copy()
    
            print(str(len(matchedDf)) + ' breakpoints were matched, ' + str(len(unmatchedDf)) + ' were not matched')
    
            matchedDf['Seq Name'] = matchedDf.index  # copy index to column 'Seq Name' in preparation for merge
            result = pd.merge(matchedDf, knownDf, on='Seq Name')  # merge matchedDf with knownDf
            result['BLAST source'] = item  # identify where the BLAST match came from in an additional column
    
            setOfFoundBPs.append(result.copy())  # save found BPs in a list of dataframes
            remainingBPs = unmatchedDf  # assign the unmatched BPs to remaining BPs
        else:
            print('no sequences remaining, BLAST search did not start')

    matchedBPDf = pd.concat(setOfFoundBPs)
    matchedBPDf.reset_index(inplace=True, drop=True) #reset index numbering after concatenation
    # write results to csv file
    matchedBPDf.to_csv(identifiedBpFilename, index=False)
    remainingBPs['Seq Name'] = remainingBPs.index  # add the breakpoint ID as a column
    remainingBPs.to_csv(unidentifiedBpFilename, index=False)
    return (matchedBPDf, remainingBPs)


def localBLASTn(querySeqListDf, subjectSeqFileName, seqCol='Unaligned', localSearch=True):
    """ 
    Do a BLASTN search given a dataframe with query sequences and fasta file with subject sequences 
       querySeqListDf is a data frame with the list of query sequences
       If using breakpoints, the dataframe should have been cleaned (cleanBreakpoints())
       seqCol is the name of the column in the querySeqListDf with the sequence strings
           CLC gives this the name 'Unaligned' by default
       subjectSeqFileName (for local BLAST) is the name of the fasta file with the subject sequences 
       subjectSeqFileName (for NCBI BLAST) is the name of the BLAST database 
       eValue is set in the global variable section
       localSearch = True when using the local BLAST search and FALSE when doing a BLAST search on the internet
    """

    # filenames
    tempBlastDb = 'tempBlastDb'  # filename and title for the temporary BLAST database,
    # the BLAST program adds its own file extensions
    querySeqName = 'tempQueryList.fa'  # filename for the temporary query sequences
    blastResultFileName = 'tempBlastResult.xml'

    # convert the input sequence list to a multi-fasta file and save it
    dataFrameToMultiFasta(querySeqListDf, seqCol=seqCol, outFileName=querySeqName)

    if localSearch:
        # ****************************
        # build local BLAST database
        # ****************************
        # first check for duplicated entries in the subject sequences.
        # This will cause the DB creation to fail otherwise
        # based on Sequence_Cleaner script (http://biopython.org/wiki/Sequence_Cleaner)
        seqIDs = []  # list of sequence IDs
        sequences = []  # list of corrected sequences
        for seq_record in SeqIO.parse(subjectSeqFileName, "fasta"):
            if seq_record.id in seqIDs:
                # print("dup found", seq_record) # for testing
                seq_record.id += "-1"  # append a -1 to duplicate records
                # print("fixed", seq_record) # for testing
            seqIDs.append(seq_record.id)
            sequences.append(seq_record)
        # subjectSeqFileName += '-1'    # for testing, to prevent overwriting the original file
        # write the corrected output to a file
        output_handle = open(subjectSeqFileName, "w")
        SeqIO.write(sequences, output_handle, "fasta")
        output_handle.close()

        # build a correctly-formatted command line request
        makeBlastDbStr = ('makeblastdb '  # command to make the database
                          '-in ' + subjectSeqFileName +  # filename of input files
                          ' -parse_seqids '  # read in sequence names
                          '-dbtype nucl '  # type of database is nucleotide
                          '-out ' + tempBlastDb +  # filename of the resulting BLAST database (6 files created)
                          ' -title ' + tempBlastDb)  # title of the resulting BLAST database
        print('\nMaking BLAST database')
        print('\t', makeBlastDbStr)  # for testing/debugging
        start = time.time()  # start a timer to keep track of progress
        run_command(makeBlastDbStr)
        # os.system(makeBlastDbStr)  # run the BLAST database creation program from the command line
        end = time.time()
        print('\t%1.4f' % (end - start), " seconds to make BLAST database")

        # *************************
        # do a local BLAST search
        # *************************
        # build a correctly-formatted command line request
        blastSearchStr = ('blastn'  # commend for the BLAST search
                          ' -out ' + blastResultFileName +  # name of the output file
                          ' -outfmt 5'  # format 5 is XML
                          ' -query ' + querySeqName +  # filename of the query sequences (multi-fasta format)
                          ' -db ' + tempBlastDb +  # filename of the BLAST database (no extensions needed)
                          ' -evalue ' + str(eValue) +  # e-value cutoff to use for the BLAST search
                          ' -num_threads 4'  # number of CPU threads to use in BLAST search
                          ' -word_size 5' # default is 11, 5 works better for short sequences
                          ' -dust no') # DUST is a program that filters low complexity sequences,
                                       # this is on by default and causes BLAST to miss some sequences
                                       # see breakpoint 2013699^2013700 near pfl in C. therm
        # run the BLAST search, use a timer to check how long it takes
        print('\nStarting BLAST search')
        start = time.time()  # start a timer to keep track of progress
        print('\t', blastSearchStr)  # for testing/debugging
        # os.system(blastSearchStr)  # run the BLAST search from the command line
        run_command(blastSearchStr)
        print('BLAST search finished')
        end = time.time()
        print('\t%1.4f' % (end - start),
              " seconds for ",
              len(querySeqListDf.index),
              " results")  # check the time to see how it's going
    else:
        # ***************************
        # do a BLAST search at NCBI
        # ***************************
        # setting megablast to True speeds things up a lot
        fasta_string = open(querySeqName).read()
        print('\nStarting BLAST search')
        start = time.time()  # start a timer to keep track of progress
        result_handle = NCBIWWW.qblast('blastn', subjectSeqFileName, fasta_string, expect=eValue, megablast=True)
        save_file = open(blastResultFileName, "w")
        save_file.write(result_handle.read())
        save_file.close()
        result_handle.close()
        end = time.time()
        print('\t%1.4f' % (end - start),
              " seconds for ",
              len(querySeqListDf.index),
              " results")  # check the time to see how it's going


    # ***********************
    # parse the results
    # ***********************
    # parse the temporary XML file
    # save the results in a dataframe
    # return the dataframe
    blast_result_table = []  # temporary list to hold 2-D table before conversion to dataframe

    print('Processing BLAST results')
    start = time.time()  # start a timer to keep track of progress
    for qresult in SearchIO.parse(blastResultFileName, 'blast-xml'):  # loop through all query sequences
        if len(qresult) > 0:  # check to see if there are any hits
            # get the first hit and first hsp (from that hit) and save results in a table
            #print('\n',seqName)
            seqName = int(qresult.id)
            queryLen = qresult.seq_len # length of the breakpoint being matched

            # loop through hits
            for hit in qresult.hits:
                seqLen = hit.seq_len
                hitID = hit.id
                first_eval = hit.hsps[0].evalue
                first_span = hit.hsps[0].query_span
                hspHitSpan = hit.hsps[0].hit_span
                hspHitStart = hit.hsps[0].hit_range[0]
                hspHitEnd = hit.hsps[0].hit_range[1]
                hspHitStrand = hit.hsps[0].hit_strand # which strand the sequence hit, either 1 or -1

                # calculate distance from start of bp
                # different calculations for Left or Right breakpoints
                # for left breakpoints, the start of the junction is at the right end of the sequence
                # for right breakpoints, the start of the junction is at the left end of the sequence
                bpName = querySeqListDf['Name'].loc[seqName]
                bpDist = 1e6 # reset value
                hitEndDist = 1e6 # reset value
                
                # calculations for left breakpoints
                if bpName == 'Left breakpoint':
                    bpDist = queryLen - hit.hsps[0].query_range[1]
                    
                    # the hit distance calculations are different for + and - strands
                    # minus strand hit distance calculations
                    if hit.hsps[0].hit_strand == -1:
                        hitEndDist = hit.hsps[0].hit_start
                    #hitEndDist = 
                    # plus strand hit distance calculations
                    else:
                        hitEndDist = (seqLen - hit.hsps[0].hit_end)
                    
                # calculations for right breakpoints
                elif bpName == 'Right breakpoint':
                    bpDist = hit.hsps[0].query_range[0]

                    # the hit distance calculations are different for + and - strands
                    # minus strand hit distance calculations
                    if hit.hsps[0].hit_strand == -1:
                        hitEndDist = (seqLen - hit.hsps[0].hit_end)
                        #if hitEndDist > 5:
                        #   print('hitEndDist = ', hitEndDist, ' strand= ', hit.hsps[0].hit_strand )
                        #  raise StopIteration
                    # plus strand hit distance calculations
                    else:
                        hitEndDist = hit.hsps[0].hit_start
                        
                # if the name isn't left or right breakpoint
                else:
                    print('------------- Error:sequence is not a breakpoint -----------------')
                    break

                # put everything together in a new row
                blast_result_row = [seqName, hitID, first_eval, bpDist, hitEndDist]
                blast_result_table.append(blast_result_row)

    blastResultDf = pd.DataFrame(blast_result_table, columns=['Seq Name',
                                                              'Hit Name',
                                                              'e-value',
                                                              'BP distance',
                                                              'Hit end distance'])
    end = time.time()
    print('\t%1.4f' % (end - start), " seconds for BLAST processing")
    return (blastResultDf)


def dataFrameToMultiFasta(inputDf, seqCol='Unaligned', outFileName = 'temp.fa'):
    """
    convert a sequence list in dataFrame format into a multi-fasta sequence file for BLAST search
    helper for local BLAST search
    """
    # dFrameIn is the dataframe
    # seqCol is the name of the column in the dataframe with the actual sequences
    # outFileName is the name of the output file

    bpSeqRecordList = [] #list to hold the sequence records before writing to file
    for index, row in inputDf.iterrows():
        tempSeqRecord = SeqRecord(Seq(row[seqCol], IUPAC.unambiguous_dna),
                                  id=str(index), #this needs to be a string for the SeqRecord function
                                  name='',
                                  description='')
        bpSeqRecordList.append(tempSeqRecord)

    # write the SeqRecord objects to a multi-fasta file in preparation for BLAST analysis
    output_handle = open(outFileName, 'w') # open a file
    SeqIO.write(bpSeqRecordList, output_handle, "fasta") # start writing breakpoints as fasta sequences
    output_handle.close() # close file at the end


def run_command(command):
    """
    Run a command on the command line and print out any messages or errors
    Sometimes the BLAST functions give error messages, and it's useful to see them
    """
    p = subprocess.Popen(command,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT)
    # print result of stdout
    if p.stdout is not None:
        print('***************START stdout***************')
        for item in p.stdout:
            print(item)
            # print(item.decode("utf-8")) #to see the output as formatted
    print('***************END stdout***************')
    # print result of stderr
    if p.stderr is not None:
        print('***************START stderr***************')
        for item in p.stderr:
            print(item)
        # print(item.decode("utf-8"))
        print('***************END stderr***************')
    return iter(p.stdout.readline, b'')
    



print('done')
