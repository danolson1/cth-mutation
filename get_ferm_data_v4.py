""" 
#####################################################################################################
This module is for getting fermentation data from the reference fermentation data excel file


Inputs(hard coded):
     - location of reference fermentation data excel file
     - value for minimum fraction of substrate consumed (good default value is 0.2)
     - list of types of media to be considered "MTC"
     
Outputs:
     - dataframe with ethanol yield data (mean, stdev, CV and count) for each strain

Dan Olson
12-12-2016
#####################################################################################################
"""

#----------------------------------------------------------------------------------------------------
# imports
#----------------------------------------------------------------------------------------------------
import os
import re
from scipy.stats import ttest_ind_from_stats
import numpy as np
import pandas as pd

#----------------------------------------------------------------------------------------------------
# constants and global variables
#----------------------------------------------------------------------------------------------------
path = r'\\jumbo.thayer.dartmouth.edu\lyndlab\biomass\C thermocellum reference fermentation data' 
fileName = r'\Compiled data through bottle 855 v4.xlsx'
filePath = path + fileName
minSubstrateFraction = 0.2 # filter out bottles that have substrate utilization lower than this
# list of 'approved' types of MTC medium (both C. therm and T. sacch
# note that I might want to get rid of the 20 g/l bottles at some point
mtcMediumList = ['MTC',
                'MTC + 0.5 g/l YE',
                'MTC + 6ug/ml thiamphenicol',
                'MTC 5 + 0.9g/L yeast extract',
                'MTC with less N2',
                'MTC+ RPMI vitamins and MEM amino acids',
                'MTC-5',
                'MTC-5 + 20g/L Cellobiose',
                'MTC-5 + 20g/L Cellobiose, 6?g/ml Thiamphenicol',
                'MTC-5 + 20g/L Cellobiose, 6?g/ml Thiamphenicol, H2 purged with H2 headspace, stationary incubation in chamber',
                'MTC-5 + 20g/L Cellobiose, 6?g/ml Thiamphenicol, N2 purged, chamber depressurized, stationary incubation in chamber',
                'MTC-5 + 20g/L Cellobiose, 6?g/ml Thiamphenicol, triple purged with chamber gas, stationary incubation in chamber',
                'MTC-5 + 20g/L Cellobiose, 6?g/ml Thiamphenicol, triple purged with chamber gas, stationary incubation in chamber, left vented in chamber',
                'MTC-5 + 6µg/ml thiamphenicol',
                'MTC-5 + 6ug/ml Thiamphenicol',
                'MTC-5, new cellobiose',
                'MTC-6',
                'MTC-6 (+4mM formate)',
                'MTC-6 low cellobiose',
                'MTC-6, 2.5 g/l cb',
                'MTC-6, 2.5 g/l cellobiose',
                'MTC-6, new cellobiose',
                'MTC-6 + 0.5g/L YE',
                'T. sacch MTC']

cleanFermData = None #global variable to hold fermentation data
rawFermData = None

                
#----------------------------------------------------------------------------------------------------
# functions
#----------------------------------------------------------------------------------------------------

def getStrainID(inStr):
    """
    get strain ID
    should start with LL or AG
    then 3-4 digits
    """
    result = ''
    p = re.compile('((AG|LL)\d{3,4})')
    if type(inStr) is str:
        regexResult = p.findall(inStr)
        if len(regexResult) > 0:
            result = regexResult[0][0]            
    return result

def setup(fp = filePath):
    """
    read excel file
    filter bottles
    get LL number
    return cleanFermData dataframe
    also save cleanFermData as a global variable """

    # read data from excel file and set up new columns
    df =pd.read_excel(fp, sheetname='bottles', skiprows=19)
    df.set_index('Bottle ID', inplace=True)
    df['substrateFracUsed'] = df['cellobiose consumed']/df['initial substrate (cellobiose)']
    df['StrainID'] = df['LL number'].apply(getStrainID) # get rid of the 'LL' from LL numbers
    global totalBottles
    totalBottles = len(df)
    
    # generate boolean filters for getting rid of rows 
    global nonRejected
    global hasEtohYield
    global hasStrainId
    global enoughSubstrateUsed
    global isMtcMedium
    
    blankRows = (df['Media'].isnull()) & (df['Blank ID'].isnull()) & (df['Bottle name'].isnull()) & (df['Strain description'].isnull())
    nonRejected = ~(df['Reject from analysis'].notnull())
    hasEtohYield = df['Yet'].notnull()
    hasStrainId = df['StrainID'] != ''
    enoughSubstrateUsed = df['substrateFracUsed'] > minSubstrateFraction
    isMtcMedium = df['Media'].isin(mtcMediumList)

    # store the raw bottle results in a dataframe for manual inspection, if needed
    global rawFermData
    rawFermData = df
    rawFermData['nonRejected'] = nonRejected
    rawFermData['hasEtohYield'] = hasEtohYield
    rawFermData['hasStrainId'] = hasStrainId
    rawFermData['enoughSubstrateUsed'] = enoughSubstrateUsed
    rawFermData['isMtcMedium'] = isMtcMedium
    
    # create a new dataframe with the non-rejected data
    df0 = df.loc[~blankRows, :]
    df1 = df0.loc[nonRejected, :]
    df2 = df1.loc[hasEtohYield, :]
    df3 = df2.loc[hasStrainId, :]
    df4 = df3.loc[enoughSubstrateUsed, :]
    df5 = df4.loc[isMtcMedium , :]

    # keep columns from 'carbon balance' to 'StrainID'
    global cleanFermData
    cleanFermData = df5.loc[:, 'carbon balance':'StrainID']
    #cleanFermData = df.loc[(nonRejected 
    #                        & hasEtohYield 
    #                        & hasStrainId 
    #                        & isMtcMedium
    #                        & enoughSubstrateUsed), 'carbon balance':'StrainID'] 

    # summarize cleaning
    print('*********** CLEANING SUMMARY ****************')
    print('initial number of bottles              = ', len(df))
    print('empty rows removed                     = ', len(df) - len(df0))
    print('number rejected                        = ', len(df0) - len(df1))
    print('rejected for no ethanol (blk bottles)  = ', len(df1) - len(df2))
    print('rejected for no strain ID              = ', len(df2) - len(df3))
    print('rejected for low substrate utilization = ', len(df3) - len(df4))
    print('rejected for non-MTC medium            = ', len(df4) - len(df5))
    print('\nfinal number of bottles                = ', len(df5))
    print('*********************************************')
    return cleanFermData


def getAgFermData(colName, fermData, parChildDf):
    """ 
    get aggregated fermentation data
    available columns include:
        carbon balance
        YpelC (yield of pellet carbon)
        Yform (yield of formate)
        Ylac (yield of lactate)
        Yac (yield of acetate)
        Yet (yield of ethanol)
        Ypyr (yield of pyruvate)
        Ymal (yield of malate)
        Yh2 (yield of hydrogen)
        PforFlux (PFOR flux = Yac + Yet - Yform)
        substrateFracUsed (fraction of initial cellobiose consumed)
    colName is the name of the column to aggreagate
    parChildDf is a dataframe with parent and child columns
        index is the strain ID as a string
        ParentID is the parent strain ID as a string"""
    
    # check to make sure colName is actually the name of a column in cleanFermData
    if colName in fermData.columns.tolist():
        aggData = fermData.groupby('StrainID').agg(['mean', 'std', 'count'])[colName]
        aggData['cv'] = aggData['std']/aggData['mean']
        # merge aggregated data with parChildDf based on parent
        result = pd.merge(parChildDf, aggData, how='left', left_on='StrainID', right_index=True)
        # merge parent strain data
    
        result = pd.merge(result, aggData, how='left', left_on='ParentID', right_index=True, suffixes=('_chi', '_par'))
        #ratioName = colName + '_ratio' # easier not to have this name change 
        ratioName = 'ratio'
        result[ratioName] = result['mean_chi'] / result['mean_par']
    
        # perform T-test based on summary statistics and save result in a new column
        # equal_var = False uses Welch's t-test
        # I think this is a correct assumption, but might want to double check
        # The resulting p-values are slightly larger than when equal_var = True
        ttestResult = pd.Series(''*len(result))
        for index, row in result.iterrows():
            #print(row)
            #print(index, ' ', row['count_chi'], ' ',row['count_par'] )
            if ((row['count_chi'] > 1) 
                 and (row['count_par'] > 1)
                 and (row['mean_chi'] > 0)
                 and (row['mean_par'] > 0)):
                (stat, pval) = ttest_ind_from_stats(row['mean_chi'], row['std_chi'], row['count_chi'], 
                                                    row['mean_par'], row['std_par'], row['count_par'],
                                                    equal_var = False)
            else:
                pval = np.nan
            ttestResult.loc[index] = pval
    
        #pvalName = colName + 'pval'
        pvalName = 'pval'
        result[pvalName] = ttestResult
        result.sort_values(pvalName, inplace=True)
    return result
        
        
