#!/usr/bin/env python3
# next remove duplicates

import pandas as pd
import pyranges as pr
import numpy as np
from utils import *

def dedup(df_all,out_dir,sample,slack=200,recipOverlap=0.8,verbose = True):
    if(verbose):
        print('deduping...',flush=True)

    # intialize variables
    df_all['id'] = df_all['sample'] + '__' + df_all['caller'] + '__' + df_all['event_id']  # uniq id to flag dups later
    df_all['isDup'] = 'N'  
    
    # generate a list of duplicated ids to mark down
    callers = sorted(df_all['caller'].unique())
    for caller in callers:
        if verbose:
            print('.' + caller + '...', end='',flush=True)
        
        df1 = df_all[df_all['caller'] == caller].copy()

        # get dups and label
        dupIDs = identifyDups(df1,slack,recipOverlap,verbose)
        df_all['isDup'] = np.where((df_all['isDup'] == 'Y') | (df_all['id'].isin(dupIDs)),'Y','N') # :|
        if verbose:
            print('done')

    # write combined with dup annotation
    outFile = out_dir + sample + '-svs.txt'
    df_all.to_csv(outFile,sep="\t",index = False)
    
    if(verbose):    
        print('deduping done')
    return(df_all)

# returns indexes of records that should be marked as duplicates
def identifyDups(df1,slack=200,recipOverlap=0.8, verbose = True):

    df1 = df1.sort_values(by='variant_id')

    df2 = convertPaired(df1,verbose)
    df2['dp'] = df2[['dp1','dp2']].min(axis=1)
    df2['spanning'] = df2[['spanning1','spanning2']].min(axis=1)            

    if verbose:
        print('searching for dups...',end='',flush=True)
    
    # ** for intrachrom events we can perform reciprical overlap **
    df2a = df2[df2['chrom1'] == df2['chrom2']].copy()
    df2a['Start'] = df2a['pos1'].astype(int) - 1; df2a['End'] = df2a['pos2'].astype(int); df2a['Chromosome'] = df2a['chrom1']

    pr2a = pr.PyRanges(df2a)
    overlap1 = pr2a.join(pr2a,report_overlap=True,how=None,preserve_order = True) # compare to self

    # infer reciprical overlap
    # the minimal ratio of overlap vs. length of compared segment - this is the max compared segment    
    df2b = overlap1.df
    df2b['width1'] = df2b['End'] - df2b['Start']
    df2b['width2'] = df2b['End_b'] - df2b['Start_b']
    df2b['maxWidth'] = df2b[['width1','width2']].max(axis=1)    
    df2b['ro'] = df2b['Overlap']/df2b['maxWidth']

    # infer differences in breakpoints
    df2b['diff1'] = abs(df2b['Start'] - df2b['Start_b'])
    df2b['diff2'] = abs(df2b['End'] - df2b['End_b'])
    df2b['maxDiff'] = df2b[['diff1','diff2']].max(axis=1)    
    
    # ** get intrachrom dups **
    df2c = df2b[(df2b['id'] != df2b['id_b']) & (df2b['ro'] >= recipOverlap) & (df2b['maxDiff'] <= slack)].copy()

    # ** add matchid so that reciprical matches so that we can prioritize and select duplicates to mark **
    df2c.loc[:,'matchID'] = df2c[['id','id_b']].apply(lambda x: '___'.join(sorted(x)),axis=1)
    df2c = df2c.sort_values(by=['dp','id'],ascending=[False,True])    # order
    idxDups = df2c['matchID'].duplicated(keep = 'first')
    df2d = df2c[idxDups]
    dupIDs1 = list(df2d['id'])  # these are intrachrom duplicates
    # **

    # ** for interchrom events we compare breakpoints separately using slack **
    df3a = df2[df2['chrom1'] != df2['chrom2']].copy()
    df3b = df3a.copy()
    df3b['Start'] = df3b['pos1'].astype(int) - 1; df3b['End'] = df3b['pos1'].astype(int); df3b['Chromosome'] = df3b['chrom1']
    df3c = df3a.copy()
    df3c['Start'] = df3c['pos2'].astype(int) - 1; df3c['End'] = df3c['pos2'].astype(int); df3c['Chromosome'] = df3b['chrom2']

    # get overlapping ranges
    pr3b = pr.PyRanges(df3b)
    pr3c = pr.PyRanges(df3c)
    overlap1 = pr3b.join(pr3b,slack = slack,report_overlap=True,how=None,preserve_order = True)
    overlap2 = pr3c.join(pr3c,slack = slack,report_overlap=True,how=None,preserve_order = True)    
    df3d = overlap1.df
    df3e = overlap2.df    

    # ** add matchid so that reciprical matches so that we can prioritize and select duplicates to mark **
    df3d1 = df3d[df3d['id'] != df3d['id_b']].copy()
    df3e1 = df3e[df3e['id'] != df3e['id_b']].copy()
    df3d1.loc[:,'matchID'] = df3d1[['id','id_b']].apply(lambda x: '___'.join(sorted(x)),axis=1)
    df3e1.loc[:,'matchID'] = df3e1[['id','id_b']].apply(lambda x: '___'.join(sorted(x)),axis=1)

    # get rid of reciprical matches
    df3d1 = df3d1.sort_values(by=['spanning','dp','id'],ascending=[False,False,True])    # order
    df3e1 = df3e1.sort_values(by=['spanning','dp','id'],ascending=[False,False,True])    # order    
    idxDups3d1 = df3d1['matchID'].duplicated(keep = 'first')
    idxDups3e1 = df3e1['matchID'].duplicated(keep = 'first')    
    df3d2 = df3d1[idxDups3d1]
    df3e2 = df3e1[idxDups3e1]        

    # get list of dups and take intersection - look for start and ends that are close via slack
    dups3d2 = list(df3d2['id'])
    dups3e2 = list(df3e2['id'])
    dupIDs2 = [id for id in dups3d2 if id in dups3e2]  # these are interchrom dups

    return(sorted(dupIDs1 + dupIDs2))
