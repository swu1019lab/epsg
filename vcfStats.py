#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import logging
import os
import sys
import time
import vcf
import pandas as pd
import numpy as np
from collections import defaultdict

__author__ = 'Xiaodong Li'
__date__   = '2021/1/25'

def vcfStats(vcfFile, outName):
    vcfReader = vcf.Reader(filename=vcfFile)
    colDict = defaultdict(list)
    for record in vcfReader:
        if not record.is_snp or len(record.alleles) != 2: sys.exit('This script only can be used for bi-allelic SNPs!')
        # For bi-allelic SNPs
        MAF = min(record.aaf[0], 1 - record.aaf[0])
        # PIC = 1 - MAF^2 - (1-MAF)^2 - 2*MAF^2*(1-MAF)^2
        PIC = 1 - np.square(MAF) - np.square(1 - MAF) - 2*np.square(MAF)*np.square(1 - MAF)
        colDict['Chr'].append(record.CHROM)
        colDict['Pos'].append(int(record.POS))
        colDict['Ref'].append(record.REF)
        colDict['Alt'].append(str(record.ALT[0]))
        colDict['Var_type'].append(record.var_type.upper())
        colDict['Var_subtype'].append(record.var_subtype)
        colDict['MAF'].append(MAF)
        colDict['PIC'].append(PIC)
        colDict['Missing'].append(1 - record.call_rate)
        colDict['Num_called'].append(record.num_called)
        colDict['Num_het'].append(record.num_het)
        colDict['Num_hom_alt'].append(record.num_hom_alt)
        colDict['Num_hom_ref'].append(record.num_hom_ref)
        colDict['Heterozygosity'].append(record.heterozygosity)
        colDict['Nucl_diversity'].append(record.nucl_diversity)

    # Construct DataFrame from dict
    df = pd.DataFrame.from_dict(colDict)
    df.to_csv("{}.csv".format(outName), index = False)
    return df
	
def main():
    parser = argparse.ArgumentParser(description='Perform basic statistics on bi-allelic SNPs from the VCF file.')
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')
    parser.add_argument("-i", "--input", type=str, required=True, help="Input VCF file")
    parser.add_argument("-o", "--output", type=str, required=True, help="Output file name prefix")
    parser.add_argument("--log", type=str, default="run.log", help="output log file")
    args = parser.parse_args()

    # Setup logging and script timing
    logging.basicConfig(filename=args.log, filemode='w', format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p', level=logging.DEBUG)
    logging.info('----- Start')
    logging.info("# {} #".format(os.path.basename(__file__)))
    timeStart = time.time()
    vcfStats(args.input, args.output)
    timeEnd = time.time()
    # End of run
    logging.info('----- Finish')
    logging.info("Elapsed time: {}".format(timeEnd - timeStart))

if __name__ == '__main__':
    # Main function
    main()
