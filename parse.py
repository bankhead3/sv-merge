#!/usr/bin/env python3
# parse vcf files

import pandas as pd
import vcf,os,re


chroms = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX', 'chrY','chrM']

# *** determine caller, parse, generate a data frame ***
def parse_vcfs(vcf_list,out_dir,sample,verbose):

    filenames = []

    # check if more than one svaba file - impacts combining calls when parsing
    multi_svaba = True if len([file for file in vcf_list if 'svaba' in file]) > 1 else False
    
    for my_vcf in vcf_list:
        # determine caller
        with open(my_vcf) as in1:
            caller = 'NA'
            for line in in1:
                if 'manta' in line or 'Manta' in line:
                    caller = 'manta'
                    break
                elif 'svaba' in line:
                    caller = 'svaba'
                    break

        # call parse function
        if caller == 'manta':
            filename1 = parse_manta(my_vcf,out_dir,sample,verbose)
            filename = modify_manta(filename1,out_dir,sample,verbose) # generated updated manta call have 2 break points per manta del,ins,dup event
        elif caller == 'svaba':
            filename = parse_svaba(my_vcf,out_dir,sample,verbose,multi_svaba)
        else:
            print('Caller not recognized!!')
            raise
        filenames.append(filename)

    filenames = sorted(list(set(filenames)))
    
    # generate combined table 
    flag = None
    for filename in filenames:
        df = pd.read_csv(filename,sep="\t",na_values = 'NA',na_filter = False)
        if flag == None:
            df_all = df.copy()
            flag = True
        else:
            df_all = pd.concat([df_all,df],axis=0)

    # return combined table of svs
    df_all = df_all.reset_index(drop = True)
    return df_all
# ***

# *** modify manta calls to be comparable with other sv callers ***
def modify_manta(inFile,out_dir,sample,verbose):
    df = pd.read_csv(inFile,sep="\t",na_values = 'NA',na_filter = False)

    # ** manta labeling patch **
    df = df.sort_values(by='variant_id')
    events = sorted(list(set(list(df['event_id']))))
    for event in events:
        idx = df['event_id'] == event
        idxs = list(df[idx].index)
        assert len(idxs) <= 2

        # for cases when we have 2 variants as part of the event
        if len(idxs) == 2:
            chrom1 = df.at[idxs[0],'chrom']
            chrom2 = df.at[idxs[1],'chrom']
            pos1 = df.at[idxs[0],'pos']
            pos2 = df.at[idxs[1],'pos']
            variant1 = df.at[idxs[0],'variant_id']
            variant2 = df.at[idxs[1],'variant_id']            

            # swap if labels if pos2 is higher - otherwise granges will complain in combine script...
            if chrom1 == chrom2 and int(pos1) > int(pos2):
                df.at[idxs[0],'variant_id'] = variant2
                df.at[idxs[1],'variant_id'] = variant1
                df.at[idxs[0],'mate_id'] = variant1
                df.at[idxs[1],'mate_id'] = variant2
    df = df.sort_values(by='variant_id')
    # ** 
    
    outFile1 = out_dir + sample + '-manta_modified.txt'
    with open(outFile1,'w') as out1:
        # assemble and write yo header
        header = ['sample','caller','event_id','variant_id','variant_type','chrom','pos','ref','alt','mate_id','tumor_discordant_rs','tumor_spanning_rs','tumor_dp','intra_chrom_event_length']
        out1.write('\t'.join(header) + '\n')

        if verbose:
            print('modifying manta calls for compatibility ...',end='')

        for index,row in df.iterrows():
            variant_type = row['variant_type']
#            print(list(row))
            if variant_type == 'BND':
                lineOut = [str(field) for field in list(row)]
                out1.write('\t'.join(lineOut) + '\n')
            elif variant_type in ['DEL','DUP','INV']:
                # write first entry
                variant_id = row['variant_id']
                row['variant_id'] = variant_id + '_bp1'
                row['mate_id'] = variant_id + '_bp2'                
                lineOut = [str(field) for field in list(row)]
                out1.write('\t'.join(lineOut) + '\n')

                # write first entry
                row['variant_id'] = variant_id + '_bp2'
                row['mate_id'] = variant_id + '_bp1'                                
                event_length = int(row['intra_chrom_event_length'])
                row['pos'] = str(int(row['pos']) + event_length)
                lineOut = [str(field) for field in list(row)]
                out1.write('\t'.join(lineOut) + '\n')
            elif variant_type == 'INS':
                # write first entry
                variant_id = row['variant_id']
                row['variant_id'] = variant_id + '_bp1'
                row['mate_id'] = variant_id + '_bp2'                
                lineOut = [str(field) for field in list(row)]
                out1.write('\t'.join(lineOut) + '\n')

                # write first entry
                row['variant_id'] = variant_id + '_bp2'
                row['mate_id'] = variant_id + '_bp1'                                
                row['pos'] = str(int(row['pos']) + 1)
                lineOut = [str(field) for field in list(row)]
                out1.write('\t'.join(lineOut) + '\n')
            else:
                print(variant_type)
                print(row)
                raise 'i no understand'

    if verbose:                
        print('done')
    return outFile1

# *** determine which of the samples (if matched) are tumor ***
# assumes higher discordant read count in the top x vcf calls
def infer_tumor_idx(vcf_reader,caller):
    num_calls = 30
    counts = {0:[],1:[]}

    call_count = 0
    for vcf_record in vcf_reader:
        if caller == 'manta':
            for idx in range(0,2):
                # first take from split read support 
                if hasattr(vcf_record.samples[idx].data,'SR'):
                    counts[idx].append(vcf_record.samples[idx].data.SR[1])
                else:
                    # else take from spanning read support
                    assert hasattr(vcf_record.samples[idx].data,'PR')
                    counts[idx].append(vcf_record.samples[idx].data.PR[1])
        elif caller == 'svaba':
            for idx in range(0,2):
                # first take from discordant read support if there is any
                if vcf_record.samples[idx].data.DR != 0:
                    counts[idx].append(vcf_record.samples[idx].data.DR)
                else:
                    # otherwise take from ad
                    counts[idx].append(vcf_record.samples[idx].data.AD)
        else:
            raise 'i no understand now to infer_tumor_idx for this caller'
        
        # read the top x calls
        call_count += 1
        if call_count >= num_calls:
            break

    # sum then and return the index that has more
    sums = {0:0,1:0}
    for my_key in sums.keys():
        for idx in range(num_calls):
            sums[my_key] += counts[my_key][idx]

    # throw a warning if second genotype is not selected
    larger_idx = 1 if sums[1] > sums[0] else 0
    if larger_idx != 1:
        print('WARNING: first vcf genotype entry is being selected as tumor due to variant read support')
        
    return larger_idx

# *** parse manta vcf and return df ***
def parse_manta(my_vcf,out_dir,sample,verbose):
    outFile1 = out_dir + sample + '-manta.txt'
    with open(outFile1,'w') as out1:
        # assemble and write yo header
        header = ['sample','caller','event_id','variant_id','variant_type','chrom','pos','ref','alt','mate_id','tumor_discordant_rs','tumor_spanning_rs','tumor_dp','intra_chrom_event_length']
        out1.write('\t'.join(header) + '\n')

        if verbose:
            print('parsing ' + my_vcf + ' ...',end='')

        # determine which genotype field to read from as tumor
        vcf_reader = vcf.Reader(filename=my_vcf)        
        assert len(vcf_reader.samples) in [1,2]
        if len(vcf_reader.samples) == 1:
            tumor_idx = 0
        else:
            # since tumor may not always be last, confirm that it is and throw an warning if it's not
            tumor_idx = infer_tumor_idx(vcf_reader, 'manta')

        # iterate and read through vcf records
        vcf_reader = vcf.Reader(filename=my_vcf)            
        for vcf_record in vcf_reader:
            record = {'sample':sample,'caller':'manta'}

            info = vcf_record.INFO

            # skip non-passing records - if PASS then filter returns empty list
            if len(vcf_record.FILTER) != 0:
                continue
            
            record['variant_id'] = vcf_record.ID            
            record['variant_type'] = info['SVTYPE']
            record['chrom'],record['pos'],record['ref'],record['alt'] = vcf_record.CHROM,vcf_record.POS,vcf_record.REF,vcf_record.ALT
            record['alt'] = str(record['alt'][0])

            # only show standard chroms
            if record['chrom'] in chroms and 'chrUn' not in record['alt'] and 'alt' not in record['alt'] and 'random' not in record['alt']:
                record['mate_id'] = 'NA' if 'MATEID' not in info else info['MATEID'][0]

                # event_id should uniquely represent each bnd pair and other events
                if 'EVENT' in vcf_record.ID:
                    record['event_id'] = vcf_record.ID
                else:
                    ids = sorted([record['variant_id'],record['mate_id']])
                    ids = [id for id in ids if id != 'NA']
                    record['event_id'] = ids[0]
                
                # get read support info for tumor only for easy support of both tumor and tumor vs normal
                if hasattr(vcf_record.samples[tumor_idx].data,'PR'):
                    values = vcf_record.samples[tumor_idx].data.PR
                    record['tumor_spanning_rs'] = values[1]
                else:
                    record['tumor_spanning_rs'] = 0
                if hasattr(vcf_record.samples[tumor_idx].data,'SR'):
                    values = vcf_record.samples[tumor_idx].data.SR
                    record['tumor_discordant_rs'] = values[1]
                else:
                    record['tumor_discordant_rs'] = 0

                # get some informatio about read depth
                # for bnd depth is provided
                if 'BND_DEPTH' in info:
                    record['tumor_dp'] = info['BND_DEPTH']
                # for everything else we use total spanning reads and discordant reads
                else:
                    values1 = vcf_record.samples[tumor_idx].data.PR
                    values2 = [0,0] if not hasattr(vcf_record.samples[tumor_idx].data,'SR') else vcf_record.samples[tumor_idx].data.SR
                    record['tumor_dp'] = int(values1[0]) + int(values1[1]) + int(values2[0]) + int(values2[1])

                # attempt to infer intra chrom length...
                record['intra_chrom_event_length'] = 'NA' if 'SVLEN' not in info else str(abs(int(info['SVLEN'][0])))
                if record['intra_chrom_event_length'] == 'NA' and record['variant_type'] == 'INS':
                    record['intra_chrom_event_length'] == 0
                elif record['intra_chrom_event_length'] == 'NA' and record['variant_type'] == 'BND':
                    result = re.findall('.*(chr[0-9XYM]+)[:]([0-9]+)',record['alt'])
                    assert len(result) == 1 and len(result[0]) == 2
                    alt_chrom,alt_pos = result[0][0],result[0][1]

                    if record['chrom'] != alt_chrom:
                        record['intra_chrom_event_length'] = -1
                    else:
                        record['intra_chrom_event_length'] = abs(int(record['pos']) - int(alt_pos))
                
                # assemble and write line out
                lineOut = []
                for field in header:
                    lineOut.append(str(record[field]))
                out1.write('\t'.join(lineOut) + '\n')
    if verbose:                
        print('done')
    return outFile1
# ***
    
# *** parse svaba sv vcf ***
def parse_svaba(my_vcf,out_dir,sample,verbose,multi_svaba):

    vcf_type = 'sv' if 'sv.vcf' in my_vcf else 'indel'
    flag = 'w' if vcf_type == 'indel' or not multi_svaba else 'a'
    
    outFile1 = out_dir + sample + '-svaba.txt'
    with open(outFile1,flag) as out1:
        header = ['sample','caller','event_id','variant_id','variant_type','chrom','pos','ref','alt','mate_id','tumor_discordant_rs','tumor_spanning_rs','tumor_dp','intra_chrom_event_length']

        # indel and sv vcfs for svaba so we need to check if file has already been created else create a header...        
        if flag == 'w':
            # assemble and write yo header
            out1.write('\t'.join(header) + '\n')
        if verbose:
            print('parsing ' + my_vcf + ' ...',end='')    
        vcf_reader = vcf.Reader(filename=my_vcf)

        # determine which genotype field to read from as tumor
        vcf_reader = vcf.Reader(filename=my_vcf)        
        assert len(vcf_reader.samples) in [1,2]
        if len(vcf_reader.samples) == 1:
            tumor_idx = 0
        else:
            # since tumor may not always be last, confirm that it is and throw an warning if it's not
            tumor_idx = infer_tumor_idx(vcf_reader, 'svaba')

        # iterate and read through vcf records
        vcf_reader = vcf.Reader(filename=my_vcf)
        for vcf_record in vcf_reader:
            record = {'sample':sample,'caller':'svaba'}

            info = vcf_record.INFO

            record['variant_type'] = 'INDEL' if vcf_type == 'indel' else info['SVTYPE']
            record['variant_id'] = vcf_record.ID
            record['chrom'],record['pos'],record['ref'],record['alt'] = vcf_record.CHROM,vcf_record.POS,vcf_record.REF,vcf_record.ALT
            record['alt'] = str(record['alt'][0])

            if record['chrom'] in chroms and 'chrUn' not in record['alt'] and 'alt' not in record['alt']:
                record['mate_id'] = 'NA' if 'MATEID' not in info else info['MATEID']

                # get event_id
                ids = sorted([record['variant_id'],record['mate_id']])
                ids = [id for id in ids if id != 'NA']
                record['event_id'] = ids[0]

                
                record['intra_chrom_event_length'] = 'NA' if 'SPAN' not in info else info['SPAN']

                # get genotype support information for tumor sample
                # leaving out normal support for now
                vcf_samples = vcf_record.samples
                vcf_sample = vcf_record.samples[tumor_idx]

                record['tumor_discordant_rs'] = vcf_sample.data.AD if not hasattr(vcf_sample.data,'DR') else vcf_sample.data.DR  # AD is used for indels where we don't have dr
                record['tumor_spanning_rs'] = vcf_sample.data.SR                
                record['tumor_dp'] = vcf_sample.data.DP

                # assemble and write line out
                lineOut = []
                for field in header:
                    lineOut.append(str(record[field]))
                out1.write('\t'.join(lineOut) + '\n')
    if verbose:                
        print('done')                
    return outFile1
# *** end svaba SV parsing ***
    
        
