#!/usr/bin/env python3
# 20231003 arb

import argparse,os
import pandas as pd
import parse as p
from parse import parse_vcfs
from compare import compare
from merge import merge
from dedup import dedup
from select import select

parser = argparse.ArgumentParser(prog='sv-merge.py', description='Combines structural variant (SV) calls from multiple caller vcfs for a given sample.', epilog='manta and svaba currently supported')
parser.add_argument('-v','--vcfs', help = 'Vcf file names separated by a comma with no spaces.  2 vcfs required for comparison. If 1 vcf is provided it will be parsed.', required = True, dest = 'vcfs')
parser.add_argument('-s','--sample-name', help = 'Sample name or label that will prefix output files',required = True, dest = 'sample')
parser.add_argument('-o','--outdir', help = 'Path to destination directory', required = False, dest = 'outdir')
parser.add_argument('--slack', help = 'Allowance in bps for imperfect variant position comparisons.  Default is 200', required = False, dest = 'slack', type = int, default = '200')
parser.add_argument('-ro','--reciprical-overlap', help = 'Proportion of overlap required for 2 intrachromosomal svs to match that is defined by the intersection divided by the larger sv length.  Default is 0 or none', required = False, dest = 'ro', type = float, default = '0')
parser.add_argument('--verbose', help = 'T or F for screen output.  Default is T for true', required = False, dest = 'verbose', default = 'T')
parser.add_argument('--caller-order', help = 'Order of variant callers to show call details when calls match. example: svaba,manta,gridss', required = False, dest = 'caller_order', default = 'svaba,manta,gridss')
parser.add_argument('-n','--number-svs', help = 'debug param for selected the first N SVs from each vcf. default is -1 to turn off and process all SVs. example: 100', required = False, dest = 'num', default = '-1')
args = parser.parse_args()

vcfs=args.vcfs
sample=args.sample
out_dir = args.outdir if args.outdir[-1] == '/' else args.outdir + '/'
slack = args.slack
verbose = True if args.verbose == 'T' else False
caller_order = args.caller_order.split(',')
numSVs = int(args.num)
recipOverlap = float(args.ro)

if not os.path.exists(out_dir):
    os.makedirs(out_dir)

def main():
    # parse vcfs
    vcf_list = sorted(vcfs.split(','))
    df_all = parse_vcfs(vcf_list,out_dir,sample,numSVs,verbose)

    # mark duplicate calls 
    df_all = dedup(df_all,out_dir,sample,slack,recipOverlap,verbose)
#    df_all = pd.read_csv('intermediate/02/results/LNCaP_APIPC-svs.txt',sep="\t")  # for testing

    if len(vcf_list) > 1:  # only possible to compare more than 1 vcf
        # compare calls
        dfMatches1 = compare(df_all,out_dir,sample,caller_order,slack,recipOverlap,verbose)

        # select matches to merge on (deal with multimatching)
#        dfMatches1 = pd.read_csv('intermediate/02/results/LNCaP_APIPC-matches1.txt',sep="\t")   # for testing
        dfMatches2 = select(dfMatches1,out_dir,sample,caller_order,verbose)

        # merge calls to be
        dfMatches2 = pd.read_csv('intermediate/02/results/LNCaP_APIPC-matches2.txt',sep="\t")   # for testing
        merge(df_all,dfMatches2,out_dir,sample,verbose,caller_order)
    
if __name__ == '__main__':
    main()
