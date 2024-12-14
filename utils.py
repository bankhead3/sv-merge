#!/usr/bin/env python3
# functions re-used in multiple modules
# note: pyranges is 0's based coords

import pandas as pd

def convertPaired(df1,verbose = True):
    if verbose:
        print('converting to paired BEs...',end='',flush = True)

    # ** create a data frame with both break points **
    id2vars = {}
    df2 = []
    for idx,row in df1.iterrows():
        lineDict = row.to_dict()

        # keep track
        id,variant_id = lineDict['id'],lineDict['variant_id']
        chrom,pos,dp,spanning = lineDict['chrom'],lineDict['pos'],lineDict['tumor_dp'],lineDict['tumor_spanning_rs']
        variantDict = {'variant_id':variant_id,'chrom':chrom,'pos':pos,'dp':dp,'spanning':spanning}
        sample,caller = lineDict['sample'],lineDict['caller']

        # add to lookup table
        if id not in id2vars:
            id2vars[id] = {}
            id2vars[id][variant_id] = variantDict
        # if see before then add to df 
        else:
            id2vars[id][variant_id] = variantDict
            variants = list(id2vars[id].keys())
            assert len(variants) == 2

            # first we order by variant id
            if variants[0] < variants[1]:
                variant1,variant2 = variants
            else:
                variant2,variant1 = variants

            # ** first we confirm that interchrom breakend order is correct by looking at chroms **
            chrom1,chrom2 = id2vars[id][variant1]['chrom'],id2vars[id][variant2]['chrom']
            if chrom1 != chrom2:
                # if interchrom sv then we order by chrom#
                chromVal1,chromVal2 = chrom1.replace('chr',''),chrom2.replace('chr','')

                # numeric chroms
                if isinstance(chromVal1,int) and isinstance(chromVal2,int):
                    chromVal1,chromVal2 = int(chromVal1),int(chromVal2)
                    if chromVal2 < chromVal1:
                        variant2,variant1 = variants
                else:
                    if chromVal2 < chromVal1:
                        variant2,variant1 = variants
                chrom1,chrom2 = id2vars[id][variant1]['chrom'],id2vars[id][variant2]['chrom']
            # **
            
            pos1,pos2 = id2vars[id][variant1]['pos'],id2vars[id][variant2]['pos']
            dp1,dp2 = id2vars[id][variant1]['dp'],id2vars[id][variant2]['dp']
            spanning1,spanning2 = id2vars[id][variant1]['spanning'],id2vars[id][variant2]['spanning']                        
            row = {'sample':sample,'caller':caller,'id':id,'chrom1':chrom1,'pos1':pos1,'dp1':dp1,'chrom2':chrom2,'pos2':pos2,'dp2':dp2,'spanning1':spanning1,'spanning2':spanning2}
            df2.append(row)
    df2 = pd.DataFrame(df2)
    return(df2)
# **
