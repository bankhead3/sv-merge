#!/usr/bin/env python3
# 20231003 arb

import argparse
import pandas as pd
import parse as p
from parse import parse_vcfs
from compare import compare
from merge import merge

parser = argparse.ArgumentParser(prog='merge-svs', description='Combines structural variant (SV) calls from multiple caller vcfs for a given sample.', epilog='manta and svaba currently supported')
parser.add_argument('-v','--vcfs', help = '2 or more vcf file names separated by a comma with no spaces', required = True, dest = 'vcfs')
parser.add_argument('-s','--sample-name', help = 'Sample name or label that will prefix output files',required = True, dest = 'sample')
parser.add_argument('-t','--save-tmp-files', help = 'Y or N indicating if intermediate temp tables should saved', required = False, dest = 'save_tmp')
parser.add_argument('-o','--outdir', help = 'Path to destination directory', required = False, dest = 'outdir')
parser.add_argument('--slack', help = 'Allowance in bps for imperfect variant position comparisons.  Default is 200', required = False, dest = 'slack', type = int)
parser.add_argument('--verbose', help = 'T or F for screen output.  Default is T', required = False, dest = 'verbose', default = 'T')
parser.add_argument('--caller-order', help = 'Order of variant callers to show call details when calls match. example: svaba,manta', required = False, dest = 'caller_order', default = 'svaba,manta')
args = parser.parse_args()

vcfs=args.vcfs
sample=args.sample
tmp_files=args.save_tmp
out_dir = args.outdir
slack = args.slack
verbose = True if args.verbose == 'T' else False
caller_order = args.caller_order.split(',')

def main():
#    print('hello world')
    vcf_list = sorted(vcfs.split(','))
    assert len(vcf_list) > 1

    # parse vcfs
    df_all = parse_vcfs(vcf_list,out_dir,sample,verbose)
    # compare calls
    df_comp = compare(df_all,out_dir,sample,slack,verbose)
    # merge calls to be 
    merge(df_all,df_comp,out_dir,sample,verbose,caller_order)
    
if __name__ == '__main__':
    main()
