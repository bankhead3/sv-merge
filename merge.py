#!/usr/bin/env python
# create final table with calls merged such that variants are unique

import numpy as np

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
            mate_id="MATEID="+str(record['mate_id'])
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
def merge(df_all,df_comp,out_dir,sample,verbose,caller_order = ['svaba','manta','gridss']):
    outFile1 = out_dir + sample + '-sv-merge.txt'
    outFile2 = out_dir + sample + '-sv-merge.vcf'    
    if verbose:
        outFiles = ','.join([outFile1,outFile2])
        print('merging and writing to ' + outFiles + ' ...',end='')    

    # filter to select SVs and matches
    df_all1 = df_all[df_all['isDup'] == 'N'].copy()    
    df_comp1 = df_comp[df_comp['isSelect'] == 'Y'].copy()

    # ** create a match lookup table **
    # ** ids in this table are representative
    # ** other ids can be sourced for each  clique
    id2clique = dict()
    clique2other = dict()    
    for idx,row in df_comp1.iterrows():
        rowDict = row.to_dict()
        id1,id2,cliqueID,other = rowDict['id'],rowDict['id_b'],rowDict['cliqueID'],rowDict['other_variant_ids']
        id2clique[id1] = cliqueID
        id2clique[id2] = cliqueID

        # include superset of matched others
        others = other.split(',')
        if cliqueID in clique2other:
            clique2other[cliqueID] = ','.join(list(sorted(set(clique2other[cliqueID].split(',') + others))))
        else:
            clique2other[cliqueID] = ','.join(others)

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

        # walk through each record and decide to write or not write
        for idx,row in df_all1.iterrows():
            record = row.to_dict()
            id = record['id']

            # ** fill in missing values **
            record['call_source'] = id.split('__')[1]
            # match
            if id in id2clique:
                record['other_variant_ids'] = clique2other[id2clique[id]]

                # fix legacy id to match
                label = record['sample'] + '__'
                record['other_variant_ids'] = record['other_variant_ids'].replace(label,'')
                
                ids1 = record['other_variant_ids'].split(',')
                callers = sorted(list(set([id.split('__')[0] for id in ids1] + [record['call_source']])))
                
                record['is_matching'] = 'Y'                
                record['num_callers'] = len(callers)
                record['callers'] = ';'.join(callers)
            # no match
            else:
                record['is_matching'] = 'N'
                record['other_variant_ids'] = record['num_callers'] = record['callers'] = 'NA'
                record['callers'] = record['call_source']
                record['num_callers'] = 1 
            # **

            if not isinstance(record['intra_chrom_event_length'],str) and not np.isnan(record['intra_chrom_event_length']):
                record['intra_chrom_event_length'] = str(int(record['intra_chrom_event_length']))
            """
            # ** testing code **
            if record['event_id'] == '1001442:1':
                print(record)
                raise
            """
            
            # if in clique but not the cliqueID then ski;
            if id in id2clique and id != id2clique[id]:
                continue
            elif id not in id2clique or (id in id2clique and id == id2clique[id]):
                # ** gather and write lineOut for txt file **
                lineOut = []
                for field in header1:
                    lineOut.append(str(record[field]))
                out1.write('\t'.join(lineOut) + '\n')

                # write vcf line out
                lineOut = get_vcf_line_out(header2,record)
                out2.write('\t'.join(lineOut) + '\n')                
