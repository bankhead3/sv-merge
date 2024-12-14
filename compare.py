#!/usr/bin/env python3
# compare normalized calls between callers and identify calls that overlap
# note: pyranges is 0's based coords

import pandas as pd
from pandas.api.types import CategoricalDtype
import pyranges as pr
from utils import *

# *** getMatches from sets of paired BEs and returns the following table ***
# sample,caller1,event_id1,chrom1,pos1,chrom2,pos2,caller2,event_id2,chrom2,pos2,diff1,diff2,ro
def getMatches(df1,df2,slack=200,recipOverlap=0.8,verbose=True):

    if verbose:
        print('matching...',end='',flush=True)
    
    # ** perform intrachrom matching using reciprical overlap and slack **
    df1a = df1[df1['chrom1'] == df1['chrom2']].copy()
    df1a['Start'] = df1a['pos1'].astype(int) - 1; df1a['End'] = df1a['pos2'].astype(int); df1a['Chromosome'] = df1a['chrom1']
    df2a = df2[df2['chrom1'] == df2['chrom2']].copy()
    df2a['Start'] = df2a['pos1'].astype(int) - 1; df2a['End'] = df2a['pos2'].astype(int); df2a['Chromosome'] = df2a['chrom1']

    # convert to rangest and compare 
    pr1a = pr.PyRanges(df1a)
    pr2a = pr.PyRanges(df2a)
    prMatch1 = pr1a.join(pr2a,report_overlap=True,how=None,preserve_order = True)  # look for any overlap

    # infer reciprical overlap
    # minimal ratio of overlap vs. length of compared segment - this is the max compared segment
    dfMatch1 = prMatch1.df
    dfMatch1['width1'] = dfMatch1['End'] - dfMatch1['Start']
    dfMatch1['width2'] = dfMatch1['End_b'] - dfMatch1['Start_b']
    dfMatch1['maxWidth'] = dfMatch1[['width1','width2']].max(axis=1)    
    dfMatch1['ro'] = dfMatch1['Overlap']/dfMatch1['maxWidth']

    # infer differences in breakpoints
    dfMatch1['diff1'] = abs(dfMatch1['Start'] - dfMatch1['Start_b'])
    dfMatch1['diff2'] = abs(dfMatch1['End'] - dfMatch1['End_b'])
    dfMatch1['maxDiff'] = dfMatch1[['diff1','diff2']].max(axis=1) # get the maximum difference between start and end breakend comparisons
    
    # apply matching criteria here to id overlaps
    dfMatch1a = dfMatch1[(dfMatch1['id'] != dfMatch1['id_b']) & (dfMatch1['ro'] >= recipOverlap) & (dfMatch1['maxDiff'] <= slack)]
    dfMatch1b = dfMatch1a[['sample','caller','id','chrom1','pos1','chrom2','pos2','caller_b','id_b','chrom1_b','pos1_b','chrom2_b','pos2_b','diff1','diff2','ro']].copy()
    # **
    
    # ** perform interchrom matching using slack only **
    # perform comparison for each partner separately
    df1b = df1[df1['chrom1'] != df1['chrom2']].copy()
    df1b['Start'] = df1b['pos1'].astype(int) - 1; df1b['End'] = df1b['pos1'].astype(int); df1b['Chromosome'] = df1b['chrom1']
    df1c = df1b.copy()
    df1c['Start'] = df1c['pos2'].astype(int) - 1; df1c['End'] = df1c['pos2'].astype(int); df1c['Chromosome'] = df1c['chrom2']

    df2b = df2[df2['chrom1'] != df2['chrom2']].copy()
    df2b['Start'] = df2b['pos1'].astype(int) - 1; df2b['End'] = df2b['pos1'].astype(int); df2b['Chromosome'] = df2b['chrom1']
    df2c = df2b.copy()
    df2c['Start'] = df2c['pos2'].astype(int) - 1; df2c['End'] = df2c['pos2'].astype(int); df2c['Chromosome'] = df2c['chrom2']

    # convert to ranges
    pr1b = pr.PyRanges(df1b)
    pr1c = pr.PyRanges(df1c)
    pr2b = pr.PyRanges(df2b)
    pr2c = pr.PyRanges(df2c)
    
    # get overlapping (matching) ranges with slack requirement
    prMatch2 = pr1b.join(pr2b,slack = slack,report_overlap=True,how=None,preserve_order = True)
    prMatch3 = pr1c.join(pr2c,slack = slack,report_overlap=True,how=None,preserve_order = True)    
    dfMatch2 = prMatch2.df
    dfMatch3 = prMatch3.df    

    # get matches in both 
    dfMatch2a = dfMatch2[dfMatch2['id'] != dfMatch2['id_b']].copy()  # remove self matches
    dfMatch3a = dfMatch3[dfMatch3['id'] != dfMatch3['id_b']].copy()  # remove self matches
    dfMatch2a.loc[:,'matchID'] = dfMatch2a[['id','id_b']].apply(lambda x: '___'.join(sorted(x)),axis=1)
    dfMatch3a.loc[:,'matchID'] = dfMatch3a[['id','id_b']].apply(lambda x: '___'.join(sorted(x)),axis=1)

    # filter to common
    matchIDs2 = set(dfMatch2a['matchID'])
    matchIDs3 = set(dfMatch3a['matchID'])
    common = matchIDs2.intersection(matchIDs3)
    dfMatch2b = dfMatch2a[dfMatch2a['matchID'].isin(common)].copy()

    # polish format
    dfMatch2b['ro'] = 'NA'
    dfMatch2b['diff1'] = abs(dfMatch2b['pos1'].astype(int) - dfMatch2b['pos1_b'].astype(int))
    dfMatch2b['diff2'] = abs(dfMatch2b['pos2'].astype(int) - dfMatch2b['pos2_b'].astype(int))
    dfMatch2c = dfMatch2b[['sample','caller','id','chrom1','pos1','chrom2','pos2','caller_b','id_b','chrom1_b','pos1_b','chrom2_b','pos2_b','diff1','diff2','ro']].copy()
    
    dfMatch = pd.concat([dfMatch1b,dfMatch2c],axis=0,ignore_index = True)
    return(dfMatch)
# ***

# *** identify overlapping svs across callers ***
def compare(df1,out_dir,sample,caller_order,slack=200,recipOverlap=0.8,verbose=True):
    callers = sorted(list(set(df1['caller'])))
    outFile1 = out_dir + sample + '-matches1.txt'    

    if verbose:
        if len(callers) == len(caller_order):
            callers = caller_order
        print('comparing ' + ' vs. '.join(callers) + ' ...', flush=True)            

    # ** first look for matching variants for calls from one caller compared to calls from other callers**
    allMatches = None
    for caller in callers:

        if verbose:
            print('.' + caller + '...',end='',flush = True)

        # filter svs to be compared for each caller
        idx_not_dups = df1['isDup'] == 'N'
        idx_caller = df1['caller'] == caller
        idx_not_caller = df1['caller'] != caller
        idx_manta_ins = (~df1['variant_id'].str.contains('INS'))  # DO NOT INCLUDE MANTA INSERTION EVENTS
        dfa1 = df1[idx_not_dups & idx_caller & idx_manta_ins].copy()        
        dfb1 = df1[idx_not_dups & idx_not_caller & idx_manta_ins].copy()                

        # ** convert to paired BEs **
        dfa2 = convertPaired(dfa1,verbose)
        dfb2 = convertPaired(dfb1,False) # set verbose to false because we already see message

        # ** match paired BEs between callers **
        matches1 = getMatches(dfa2,dfb2,slack,recipOverlap,verbose)

        # concat matches from other callers
        allMatches = matches1 if not isinstance(allMatches,pd.DataFrame) else pd.concat([allMatches,matches1],axis=0,ignore_index = True)
        
        if verbose:
            print('done')

    # ** for scenario when no overlapping variant_ids have been found
    if not isinstance(allMatches,pd.DataFrame):
        lineOut = 'no matching variants identified!'
        with open(outFile1,'w') as out1:
            out1.write(lineOut + '\n')
        dfMatches = pd.DataFrame({'event_id':[lineOut]})
        if verbose:
            print('NO MATCHES FOUND!...done')
        return dfMatches
    # **

    # ** remove reciprical matches selecting by caller order ** 
    allMatches.loc[:,'matchID'] = allMatches[['id','id_b']].apply(lambda x: '___'.join(sorted(x)),axis=1)
    callCategory = CategoricalDtype(categories = caller_order,ordered = True)
    allMatches['caller'] = allMatches['caller'].astype(callCategory)
    allMatches = allMatches.sort_values(by='caller')
    idxDups = allMatches['matchID'].duplicated(keep = 'last')
    matches2 = allMatches[idxDups]

    matches2.to_csv(outFile1,sep="\t",index = False)

    if verbose:
        print('done')
        
    return matches2
