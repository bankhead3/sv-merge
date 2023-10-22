#!/usr/bin/env python
# create final table with calls merged such that variants are unique

def merge(df_all,df_comp,out_dir,sample,verbose,caller_order = ['svaba','manta']):
    outFile1 = out_dir + sample + '-sv-merge.txt'
    if verbose:
        print('merging and writing to ' + outFile1 + '...',end='')    

    # get list of variant universe
    df_all['unique_variant_id'] = df_all['caller'] + '__' + df_all['variant_id']
    df_comp['unique_variant_id_matches'] = df_comp['caller2'] + '__' + df_comp['variant_id2']        
    unique_variant_ids = sorted(list(df_all['unique_variant_id']))

    written_variant_ids = dict() # keep track of variants written so we don't write twice
    
    with open(outFile1,'w') as out1:

        # assemble and write header
        header = ['sample','event_id','is_matching','variant_id','chrom','pos','ref','alt','mate_id','callers','num_callers','call_source','tumor_discordant_rs','tumor_spanning_rs','tumor_dp','intra_chrom_event_length','other_variant_ids']
        out1.write('\t'.join(header) + '\n')
        
        # get variant universe
        for unique_variant_id in unique_variant_ids:
            caller,variant_id = unique_variant_id.split('__')
            
            # check that it hasn't already been written
            if unique_variant_id in written_variant_ids:
                continue

            # look up event
            idx_all = (df_all['variant_id'] == variant_id) & (df_all['caller'] == caller)
            assert sum(idx_all) == 1
            event_id = list(df_all[idx_all]['event_id'])
            assert len(event_id) == 1

            # find out if there is a match in comp table
            idx_comp_variant = (df_comp['variant_id'] == variant_id) & (df_comp['is_event_match'] == True) 

            # ** category #1 - we have a match **
            if sum(idx_comp_variant) >= 1:

                # * variant_id - pull values from priortized variant caller *
                # identify matching records to select from across callers
                unique_variant_id_matches = list(df_comp[idx_comp_variant]['unique_variant_id_matches'])
                unique_variant_id_candidates = [unique_variant_id] + unique_variant_id_matches

                # identify which one should be sourced from
                select_uvid = 'NA'                
                for caller_ in caller_order:
                    select_uvids = sorted([id for id in unique_variant_id_candidates if caller_ in id and id not in written_variant_ids])
                    if len(select_uvids) >= 1:
                        select_uvid = select_uvids[0]
                        break
                """
                # for testing...
                if variant_id in ['871774:1','871803:1','MantaBND:113677:0:2:0:0:0:1']:
                    print(variant_id)
                    print(unique_variant_id_matches)
                    print(unique_variant_id_candidates)
                    print(select_uvid)
                    raise
                """
                    
                record = df_all[df_all['unique_variant_id'] == select_uvid].squeeze()
                record = record.to_dict()

                # update with missing columns
                record['is_matching'] = 'Y'                
                record['callers'] = ';'.join(sorted(list(set([uvid.split('__')[0] for uvid in unique_variant_id_candidates]))))
                record['num_callers'] = len(unique_variant_id_candidates)
                record['call_source'] = select_uvid.split('__')[0]
                record['other_variant_ids'] = ';'.join([id for id in unique_variant_id_candidates if id != select_uvid])                                    

                # gather and write lineOut
                lineOut = []
                for field in header:
                    lineOut.append(str(record[field]))
                out1.write('\t'.join(lineOut) + '\n')

                # keep track of variants written
                for uvid in unique_variant_id_candidates:
                    written_variant_ids[uvid] = '.'
                    
            # ** category #2 - no match **
            else:
                # * variant_id - pull values from df_all from original caller * 
                record = df_all[idx_all].squeeze()
                record = record.to_dict()

                # update with missing columns
                record['is_matching'] = 'N'                
                record['callers'] = caller
                record['num_callers'] = 1
                record['call_source'] = caller
                record['other_variant_ids'] = '.'
                
                # gather and write lineOut
                lineOut = []
                for field in header:
                    lineOut.append(str(record[field]))
                out1.write('\t'.join(lineOut) + '\n')

                # keep track of variants written
                written_variant_ids[unique_variant_id] = '.'
                
    if verbose:
        print('done')

    
