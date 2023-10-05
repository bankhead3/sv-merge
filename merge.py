#!/usr/bin/env python
# create final table with calls merged such that variants are unique


def merge(df_all,df_comp,out_dir,sample,verbose,caller_order = ['svaba','manta']):
    outFile1 = out_dir + sample + '-sv-merge.txt'
    if verbose:
        print('merging and writing to ' + outFile1 + '...',end='')    

    # get list of variant universe
    df_all['unique_variant_id'] = df_all['caller'] + '__' + df_all['variant_id']
    df_comp['unique_variant_id_matches'] = df_comp['caller2'] + '__' + df_comp['variant_id2']        
    unique_variant_ids = df_all['unique_variant_id']
    
    written_variant_ids = dict()
    
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

            # look up event and mate
            idx_all = (df_all['variant_id'] == variant_id) & (df_all['caller'] == caller)
            assert sum(idx_all) == 1          
            event_id = list(df_all[idx_all]['event_id'])
            assert len(event_id) == 1
            idx_event = df_all['event_id'] == event_id[0]
            assert sum(idx_event) in [1,2]
            mate_id = list(df_all[idx_all]['mate_id'])[0]
        
            # find out if there is a match in comp table
            idx_comp_variant = (df_comp['variant_id'] == variant_id) & (df_comp['is_partial_event_match'] == True)
            idx_comp_mate = (df_comp['variant_id'] == mate_id) & (df_comp['is_partial_event_match'] == True)        

            # ** category #1 - we have a match **
            if sum(idx_comp_variant) >= 1 and (sum(idx_comp_mate) >= 1 or mate_id == 'NA'):
                # * variant_id - pull values from priortized variant caller *
                # identify matching records to select from 
                unique_variant_id_matches = list(df_comp[idx_comp_variant]['unique_variant_id_matches'])

                unique_variant_id_candidates = [unique_variant_id] + unique_variant_id_matches
                select_uvid = select_caller = 'NA'
                for caller_ in caller_order:
                    select_uvids = [id for id in unique_variant_id_candidates if caller_ in id and id not in written_variant_ids]
                    if len(select_uvids) == 1:
                        select_uvid = select_uvids[0]
                        select_caller = caller_
                        break

                record = df_all[df_all['unique_variant_id'] == select_uvid].squeeze()
                record = record.to_dict()

                # update with missing columns
                record['is_matching'] = 'Y'                
                record['callers'] = ';'.join(sorted([uvid.split('__')[0] for uvid in unique_variant_id_candidates]))
                record['num_callers'] = len(unique_variant_id_candidates)
                record['call_source'] = select_uvid.split('__')[0]
                record['other_variant_ids'] = ';'.join([id for id in unique_variant_id_candidates if id != select_uvid])                                    

                # gather and write lineOut
                lineOut = []
                for field in header:
                    lineOut.append(str(record[field]))
                out1.write('\t'.join(lineOut) + '\n')

                # * mate_id - pull values from priortized variant caller *
                # identify matching records to select from  # assume that if mate is NA then this holds for matches
                unique_mate_id_candidates = []
                if mate_id != 'NA' and caller + '__' + mate_id not in written_variant_ids:
                    unique_mate_id_matches = list(df_comp[idx_comp_mate]['unique_variant_id_matches'])                

                    unique_mate_id_candidates = [caller + '__' + mate_id] + unique_mate_id_matches
                    select_umid = 'NA'

                    # using same caller as for variant_id
                    select_umids = [id for id in unique_mate_id_candidates if select_caller in id and id not in written_variant_ids]
                    # it may be we already gave the mate away for a given caller
                    if len(select_umids) == 0:
                        for caller_ in caller_order:
                            select_umids = [id for id in unique_mate_id_candidates if caller_ in id and id not in written_variant_ids]                        
                            if len(select_umids) == 1:
                                break                    
                    assert len(select_umids) == 1

                    select_umid = select_umids[0]                    
                    record = df_all[df_all['unique_variant_id'] == select_umid].squeeze()
                    record = record.to_dict()

                    # update with missing columns
                    record['is_matching'] = 'Y'
                    record['callers'] = ';'.join(sorted([umid.split('__')[0] for umid in unique_mate_id_candidates]))
                    record['num_callers'] = len(unique_mate_id_candidates)
                    record['call_source'] = select_umid.split('__')[0]
                    record['other_variant_ids'] = ';'.join([id for id in unique_mate_id_candidates if id != select_umid])                    

                    # gather and write lineOut
                    lineOut = []
                    for field in header:
                        lineOut.append(str(record[field]))
                    out1.write('\t'.join(lineOut) + '\n')

                # record matching variant and mate_id of matching callers
                combined_uvids = unique_variant_id_candidates + unique_mate_id_candidates
                for uvid in combined_uvids:
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

                # * mate_id - pull values from df_all from original caller *                 
                if mate_id != 'NA' and caller + '__' + mate_id not in written_variant_ids:
                    idx_all = (df_all['variant_id'] == mate_id) & (df_all['caller'] == caller)

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

            written_variant_ids[unique_variant_id] = '.'
            written_variant_ids[caller + '__' + mate_id] = '.'
    if verbose:
        print('done')

    
