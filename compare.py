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
        df_joined['is_variant_match'] = True

        if first:
            df_comp = df_joined
            first = False
        else:
            df_comp = pd.concat([df_comp,df_joined],axis=0)            

    # ** determine if mate ids match **
    df_comp = df_comp.reset_index(drop = True)
    df_comp['is_event_match'] = False
    for index,row in df_comp.iterrows():
        variant_id = row['variant_id']
        mate_id = row['mate_id']
        mate_id2 = row['mate_id2']        
        event_id = row['event_id']
        event_id2 = row['event_id2']        
        is_variant_match = row['is_variant_match']
        variant_type = row['variant_type']
        variant_type2 = row['variant_type2']        

        # check number of events associated with variant
        idx_event_matches = (df_comp['event_id'] == event_id) & (df_comp['is_variant_match'] == True) & (df_comp['event_id2'] == event_id2)
        idx_variant = (df_comp['variant_id'] == variant_id) & (df_comp['is_variant_match'] == True) & idx_event_matches
        idx_mate = (df_comp['variant_id'] == mate_id) & (df_comp['is_variant_match'] == True) & idx_event_matches
        num_event_variant_matches = sum(idx_event_matches)

        # need to make sure that variant and mate map to different variant_ids
        matching_variant_id = list(set(df_comp['variant_id2'][idx_variant]))
        matching_variant_id2 = list(set(df_comp['variant_id2'][idx_mate]))
        ambig = True if len(matching_variant_id) == 1 and len(matching_variant_id2) == 1 and matching_variant_id[0] == matching_variant_id2[0] else False

        # strict bnd to bnd comparison does not work well between manta and svaba since manta represents many svs as single bps in vcf
        if is_variant_match and num_event_variant_matches == 2 and sum(idx_variant) >= 1 and sum(idx_mate) >= 1 and not ambig:
            df_comp.loc[index,'is_event_match'] = True
        """
        # for troubleshooting
        if event_id == "1043085:1":
            print(variant_id)
            print(num_event_variant_matches)
            print(sum(idx_variant))
            print(sum(idx_mate))
            print(df_comp.loc[index,])
            print('.')
        """

    outFile1 = out_dir + sample + '-comparison.txt'    
    df_comp.to_csv(outFile1,sep="\t",index = False)
    if verbose:
        print('done')
    
    return df_comp

