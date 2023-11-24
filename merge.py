#!/usr/bin/env python
# create final table with calls merged such that variants are unique

# ** create vcf formatted line out **
def get_vcf_line_out(header2,record):
    lookup = {'#CHROM':'chrom','POS':'pos','REF':'ref','ALT':'alt','ID':'variant_id','QUAL':'100','FILTER':'PASS'}
    lineOut = []    
    for field in header2:
        if field in lookup and lookup[field] in record:
            lineOut.append(str(record[lookup[field]]))
        elif field in lookup:
            lineOut.append(str(lookup[field]))
        elif field == 'INFO':
            span = "SPAN="+str(record['intra_chrom_event_length'])
            mate_id="MATEID="+record['mate_id']
            id_="ID="+record['variant_id']
            event="event="+record['event_id']
            info = ';'.join([span,mate_id,id_,event])
            lineOut.append(info)
        elif field == 'FORMAT':
            format_ = ':'.join(["GT","AD","DP","DR","SR"])
            lineOut.append(format_)
        elif field == 'tumor_sample':
            disc = record['tumor_discordant_rs']
            span = record['tumor_spanning_rs']
            ad = disc if int(disc) > int(span) else span
            tumor_sample = ':'.join(['0/1',str(ad),str(record['tumor_dp']),str(disc),str(span)])
            assert format_.count(':') == tumor_sample.count(':')
            lineOut.append(tumor_sample)
        else:
            print('i no understand')
            print(field)
            raise

    return lineOut
# **

# ** combine together matching calls and write to a sensible tab-delimited file and a non-sensible vcf file **
def merge(df_all,df_comp,out_dir,sample,verbose,caller_order = ['svaba','manta']):
    outFile1 = out_dir + sample + '-sv-merge.txt'
    outFile2 = out_dir + sample + '-sv-merge.vcf'    
    if verbose:
        outFiles = ','.join([outFile1,outFile2])
        print('merging and writing to ' + outFiles + ' ...',end='')    

    # get list of variant universe
    df_all['unique_variant_id'] = df_all['caller'] + '__' + df_all['variant_id']
    df_all['unique_event_id'] = df_all['caller'] + '__' + df_all['event_id']    
    df_comp['unique_variant_id_matches'] = df_comp['caller2'] + '__' + df_comp['variant_id2']

    # re-index for faster lookup
    df_all_idx = df_all.set_index('unique_variant_id')
    df_all_idx.sort_index(inplace = True)
    
    written_variant_ids = dict() # keep track of variants written so we don't write twice
    
    with open(outFile1,'w') as out1, open(outFile2,'w') as out2:

        # assemble and write header
        header1 = ['sample','event_id','is_matching','variant_id','chrom','pos','ref','alt','mate_id','callers','num_callers','call_source','tumor_discordant_rs','tumor_spanning_rs','tumor_dp','intra_chrom_event_length','other_variant_ids']
        out1.write('\t'.join(header1) + '\n')

        # add vcf header
        header_rows = ["##fileformat=VCFv4.2",
                       "##filedate=20231122",
                       '##reference=/fh/fast/ha_g/grp/reference/GRCh38/Broad_bundle_hg38/v0/Homo_sapiens_assembly38.fasta',
                       '##INFO=<ID=MATEID,Number=1,Type=String,Description="ID of mate breakends">',
                       '##INFO=<ID=ID,Number=1,Type=String,Description="variant_id">',
                       '##INFO=<ID=event,Number=1,Type=String,Description="event_id">',                                              
                       '##INFO=<ID=SPAN,Number=1,Type=Integer,Description="Distance between the breakpoints. -1 for interchromosomal">',
                       '##FORMAT=<ID=AD,Number=1,Type=Integer,Description="Allele depth: Number of reads supporting the variant">',
                       '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth of coverage: Number of reads covering site.">',
                       '##FORMAT=<ID=GT,Number=1,Type=String,Description="Most likely genotype">',
                       '##FORMAT=<ID=DR,Number=1,Type=Integer,Description="Number of discordant-supported reads for this variant">',
                       '##FORMAT=<ID=SR,Number=1,Type=Integer,Description="Number of spanning reads for this variants">',
                       "##SAMPLE=<ID=tumor_sample>"]
        for row in header_rows:
            out2.write(row + '\n')                
        header2 = ['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','tumor_sample']
        out2.write('\t'.join(header2) + '\n')        

        # walk through event variant ids and lookup matches in comp table
#        count = 0
        grouped_df = df_all.groupby('unique_event_id')
        for event_id,group in grouped_df:
            for idx,row in group.iterrows():
#                count += 1  # temp
#                print(count) # temp

                unique_variant_id = row['unique_variant_id']
                caller,variant_id = unique_variant_id.split('__')

                # check that it hasn't already been written and that it has orientation information
                if unique_variant_id in written_variant_ids or not ('[' in row['alt'] or ']' in row['alt']):
                    continue

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

                    record = df_all_idx.loc[[select_uvid],:].squeeze()  # use indexing for faster lookup
                    record = record.to_dict()

                    # update with missing columns
                    record['is_matching'] = 'Y'                
                    record['callers'] = ';'.join(sorted(list(set([uvid.split('__')[0] for uvid in unique_variant_id_candidates]))))
                    record['num_callers'] = len(unique_variant_id_candidates)
                    record['call_source'] = select_uvid.split('__')[0]
                    record['other_variant_ids'] = ';'.join([id for id in unique_variant_id_candidates if id != select_uvid])                                    

                    # gather and write lineOut
                    lineOut = []
                    for field in header1:
                        lineOut.append(str(record[field]))
                    out1.write('\t'.join(lineOut) + '\n')

                    # write vcf line out
                    lineOut = get_vcf_line_out(header2,record)
                    out2.write('\t'.join(lineOut) + '\n')                
                
                    # keep track of variants written
                    for uvid in unique_variant_id_candidates:
                        written_variant_ids[uvid] = '.'
                    
                # ** category #2 - no match **
                else:
                    # * variant_id - pull values from df_all from original caller * 
                    record = row
                    record = record.to_dict()

                    # update with missing columns
                    record['is_matching'] = 'N'                
                    record['callers'] = caller
                    record['num_callers'] = 1
                    record['call_source'] = caller
                    record['other_variant_ids'] = '.'
                
                    # gather and write lineOut
                    lineOut = []
                    for field in header1:
                        lineOut.append(str(record[field]))
                    out1.write('\t'.join(lineOut) + '\n')

                    # write vcf line out
                    lineOut = get_vcf_line_out(header2,record)
                    out2.write('\t'.join(lineOut) + '\n')                
                
                    # keep track of variants written
                    written_variant_ids[unique_variant_id] = '.'

    if verbose:
        print('done')
# **
    
