#!/usr/bin/env python3
# compare normalized calls between callers and identify calls that overlap
# note: pyranges is 0's based coords

import pandas as pd
import pyranges as pr

# *** identify overlapping svs across callers ***
def compare(df1,out_dir,sample,slack=200,verbose=True):
    callers = sorted(list(set(df1['caller'])))
    
    if verbose:
        print('comparing ' + ' vs. '.join(callers) + ' ...',end='')

    # ** first look for matching variants **
    first = True
    for caller in callers:

        idx_caller = df1['caller'] == caller
        idx_not_caller = df1['caller'] != caller
        idx_manta_ins = (~df1['variant_id'].str.contains('INS'))  # DO NOT INCLUDE MANTA INSERTION EVENTS
        df1a = df1[['sample','chrom','pos','variant_id','event_id','mate_id','caller','variant_type']][idx_caller & idx_manta_ins].copy()
        df1b = df1[['chrom','pos','variant_id','event_id','mate_id','caller','variant_type']][idx_not_caller & idx_manta_ins].copy()        

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

    
    # ** next look for matching events consisting of matching variant pairs **
    df_comp = df_comp.reset_index(drop = True)
    df_comp['is_event_match'] = False
    event_ids = sorted(list(df_comp['event_id']))
    for event_id in event_ids:
#        print(event_id)

        idx_event = df_comp['event_id'] == event_id        
        variants = sorted(set(list(df_comp['variant_id'][idx_event])))

        # require: both variants of an event match 
        if len(variants) != 2:
            continue

        variant1,variant2 = variants
        idx_variant1 = df_comp['variant_id'] == variant1
        idx_variant2 = df_comp['variant_id'] == variant2
        matching_events = sorted(list(set(df_comp['event_id2'][idx_variant1 | idx_variant2])))
        matching_variants = sorted(list(set(df_comp['variant_id2'][idx_variant1 | idx_variant2])))        

        """
        # uncomment to qc
        if event_id == '799314:1':
            print('QC HERE')
            print(variants)
            print(matching_events)
        """
        # matching sceniario #0: both variants map to the same matching variant
        if len(matching_events) == 1 and len(variants) == 2 and len(matching_variants) == 1:
            continue
        # matching scenario #1: perfect event matching were both event variants match both varaints of another caller event
        elif len(matching_events) == 1 and len(variants) == 2 and len(matching_variants) == 2 and sum(idx_variant1 | idx_variant2) == 2:
            df_comp.loc[idx_variant1,'is_event_match'] = True
            df_comp.loc[idx_variant2,'is_event_match'] = True            
        elif len(matching_events) > 1 and len(variants) == 2:
            # matching scenario #2: variants of an event match multiple events
            for matching_event in matching_events:
#                print('.' + matching_event)

                # get matching variants 
                idx_variant1_match = (df_comp['event_id2'] == matching_event) & idx_variant1
                idx_variant2_match = (df_comp['event_id2'] == matching_event) & idx_variant2               
                variant1_matches = list(df_comp['variant_id2'][idx_variant1_match])
                variant2_matches = list(df_comp['variant_id2'][idx_variant2_match])

                # scenario 2a: best case scenario both variants match only once to separate events
                if len(variant1_matches) == 1 and len(variant2_matches) == 1 and variant1_matches[0] != variant2_matches[0]:
                    df_comp.loc[idx_variant1_match,'is_event_match'] = True
                    df_comp.loc[idx_variant2_match,'is_event_match'] = True
                # scenario 2b: event variants match the same event
                elif len(variant1_matches) == 1 and len(variant2_matches) == 1 and variant1_matches[0] == variant2_matches[0]:
                    continue
                # scenario 2c: if only 1 matching variant then is not a matching event                                    
                elif len(variant1_matches) == 0 or len(variant2_matches) == 0:
                    continue
                # scenario 2d: multimatching in the same event (nasty!)
                elif len(variant1_matches) > 1 or len(variant2_matches) > 1:
#                    print('nasty')                    
                    # select the closest matches
                    differences1 = list(df_comp['difference'][idx_variant1_match])
                    differences2 = list(df_comp['difference'][idx_variant2_match])
                    min_difference = min(differences1 + differences2)

                    if min_difference in differences1:
                        idx = differences1.index(min_difference)
                        variant1_match = variant1_matches[idx]
                        variant2_match = [variant for variant in variant2_matches if variant != variant1_match][0]
                    else:
                        idx = differences2.index(min_difference)
                        variant2_match = variant2_matches[idx]
                        variant1_match = [variant for variant in variant1_matches if variant != variant2_match][0]

                    # set only the best matching pair as a matching event
                    idx_variant1_match_unique = idx_variant1_match & (df_comp['variant_id2'] == variant1_match)
                    idx_variant2_match_unique = idx_variant2_match & (df_comp['variant_id2'] == variant2_match)
                    assert sum(idx_variant1_match_unique) == 1 and sum(idx_variant2_match_unique) == 1
                    df_comp.loc[idx_variant1_match_unique,'is_event_match'] = True
                    df_comp.loc[idx_variant2_match_unique,'is_event_match'] = True
                else:
                    print(variant1_match)
                    print(variant2_match)                
                    raise 'unexpected event matching scenario'

    outFile1 = out_dir + sample + '-comparison.txt'    
    df_comp.to_csv(outFile1,sep="\t",index = False)
    if verbose:
        print('done')
    
    return df_comp
