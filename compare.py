#!/usr/bin/env python3
# compare normalized calls between callers and identify calls that overlap
# selects allos for slack when comparing
# note: pyranges is 0's based coords

import pandas as pd
import pyranges as pr

# *** identify overlapping svs across callers ***
def compare(df1,out_dir,sample,slack=200,verbose=True):
    callers = sorted(list(set(df1['caller'])))
    
    if verbose:
        print('comparing ' + ' vs. '.join(callers) + ' ...',end='')

    # ** first look for matching variants and select the best match **
    first = True
    for caller in callers:

        idx_caller = df1['caller'] == caller
        idx_not_caller = df1['caller'] != caller        
        df1a = df1[['sample','chrom','pos','variant_id','event_id','mate_id','caller','variant_type']][idx_caller].copy()
        df1b = df1[['chrom','pos','variant_id','event_id','mate_id','caller','variant_type']][idx_not_caller].copy()        


        df1a['Start'] = df1a['pos'] - 1; df1a['End'] = df1a['pos']; df1a['Chromosome'] = df1a['chrom']
        df1b['Start'] = df1b['pos'] - 1; df1b['End'] = df1b['pos']; df1b['Chromosome'] = df1b['chrom']        
        
        pr1a = pr.PyRanges(df1a)
        pr1b = pr.PyRanges(df1b)        
        
        pr_joined = pr1a.join(pr1b,slack=slack,how = None, suffix='2',preserve_order=True)

        # add difference
        df_joined = pr_joined.df
        df_joined['difference'] = abs(df_joined['pos'] - df_joined['pos2'])
        
        # polish columns 
        dfj_sorted = df_joined.sort_values(by=['caller','difference','pos2','variant_id2'], ascending = [True,True,True,True])

        idx_first_duplicate = ~ dfj_sorted['variant_id'].duplicated(keep='first')
        dfj_sorted['is_variant_match'] = idx_first_duplicate

        dfjs_polished = dfj_sorted.sort_values(by=['chrom','pos','difference'], ascending = [True,True,True])
        dfjs_polished.drop(['Chromosome','Start','End','Start2','End2'],inplace=True,axis=1)
        dfjs_polished.reset_index(drop = True)

        if first:
            df_comp = dfjs_polished
            first = False
        else:
            df_comp = pd.concat([df_comp,dfjs_polished],axis=0)

    # ** determine if mate ids match **
    df_comp = df_comp.reset_index(drop = True)
    df_comp['is_event_match'] = False
    df_comp['is_partial_event_match'] = False    
    for index,row in df_comp.iterrows():
        variant_id = row['variant_id']
        mate_id = row['mate_id']
        mate_id2 = row['mate_id2']        
        event_id = row['event_id']
        is_variant_match = row['is_variant_match']
        variant_type = row['variant_type']
        variant_type2 = row['variant_type2']        

        # check number of events associated with variant 
        idx_event_matches = (df_comp['event_id'] == event_id) & (df_comp['is_variant_match'] == True)
        num_event_variant_matches = sum(idx_event_matches)
        assert num_event_variant_matches in [1,2]
        num_unique_event_matches = len(list(set(df_comp[idx_event_matches]['event_id2'])))                

        # strict bnd to bnd comparison does not work well between manta and svaba since manta represents many svs as single bps in vcf
        if is_variant_match and mate_id == 'NA' and mate_id2 == 'NA' and num_event_variant_matches == 1:
            df_comp.loc[index,'is_event_match'] = True            
        elif is_variant_match and mate_id != 'NA' and mate_id2 != 'NA' and num_event_variant_matches == 2 and num_unique_event_matches == 1:
            df_comp.loc[index,'is_event_match'] = True

        # less stringent comparison just requires that variants in the vcf match up for a given event 
        if is_variant_match and mate_id == 'NA':
            df_comp.loc[index,'is_partial_event_match'] = True
        elif is_variant_match and mate_id != 'NA' and num_event_variant_matches == 2 and num_unique_event_matches >= 1:
            df_comp.loc[index,'is_partial_event_match'] = True

    outFile1 = out_dir + sample + '-comparison.txt'    
    df_comp.to_csv(outFile1,sep="\t",index = False)
    if verbose:
        print('done')
    
    return df_comp

