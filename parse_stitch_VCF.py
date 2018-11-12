#!/usr/bin/env python

import os
import sys
import vcf
import pprint
import logging
import argparse
import warnings
from collections import OrderedDict


def openIOFile(fileName, path = '', mode = 'r'):
    filePath = os.path.join(path, fileName)
    try:
        logging.debug(f'trying to open file [{filePath}]')
        openedFile = open(filePath, mode)
        logging.debug(f'opened file [{filePath}]')
    except IOError:
        logging.critical(f'The file [{filePath}] isn\'t accessible.')
        sys.exit(0)

    return openedFile

def verifyParentalGenosFileStructure(headerLine):
    if headerLine[:4] == ['CHROM', 'POS', 'REF', 'ALT']:
        result = True
    else:
        result = False

    return result

def extractParentGenosForGivenChrom(openedFile, chrom, header):
    genoDict = OrderedDict()
    for line in openedFile:
        elements = line.strip().split()
        logging.info(elements)
        if elements[0] == chrom:
            tempDict = OrderedDict(zip(header, elements))
            logging.info(f'\n{pprint.pformat(tempDict)}\n')

            genoDict[tempDict['POS']] = {k:tempDict[k] for k in header[2:]}
    logging.debug(f'\n{pprint.pformat(genoDict)}\n')

    try:
        assert len(genoDict.items()) > 0
    except AssertionError:
        message = f'There are no entries for chromosome {chrom} in the parental genotypes file.'
        logging.critical(message)
        sys.exit(message)

    return genoDict

def main():
    """

    :return:
    """

    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-v', '--inVCF',
        required=True,
        help='The input VCF file name including full/relative path')
    parser.add_argument(
        '-f', '--founders',
        required=True,
        help='The parental genotypes in TSV format [CHROM, POS, REF, ALT, GT1, GT2]')
    parser.add_argument(
        '-c', '--chrom',
        required=True,
        help='specify which chromosome the VCF file is from')
    parser.add_argument(
        '-p', '--prefix',
        default='out',
        help='the prefix for the output files')
    parser.add_argument(
        '-o', '--outputDir',
        default='',
        help='the name of the output directory')
    args = parser.parse_args()

    logger = logging.getLogger('root')
    FORMAT = "[%(filename)s:%(lineno)s - %(funcName)20s() ] %(levelname)s - %(message)s"
    logging.basicConfig(level=logging.INFO, format=FORMAT)

    parentalGenosDict = OrderedDict()
    with openIOFile(args.founders) as parentalGenosInput:
        parentalGenosHeader = parentalGenosInput.readline().strip().split()
        logging.info(f'header: {parentalGenosHeader}')

        try:
            assert verifyParentalGenosFileStructure(parentalGenosHeader)
        except AssertionError:
            message = f'the header of the parental genotypes file [{args.founders}] is missing. ' \
                      f'Add a header: CHROM, POS, REF, ALT, Parent1_Geno, Parent2_Geno'
            logging.critical(message)
            sys.exit(message)

        parents = parentalGenosHeader[4:6]
        logging.info(f'parents: {parents}')
        parentalGenosDict = extractParentGenosForGivenChrom(parentalGenosInput, args.chrom, parentalGenosHeader)



    vcfFileInput = openIOFile(args.inVCF)
    vcf_reader = vcf.Reader(vcfFileInput)

    sample_names = vcf_reader.samples

    # outputFN = f'{args.prefix}.{args.chrom}.genos.tsv'
    # outputFile = openIOFile(outputFN, args.outputDir, 'w')
    # outputLine = '\t'.join(map(str, header))


    for record in vcf_reader:
        # print(record.CHROM, record.POS, record.num_called, record.num_unknown)

        # try:
        #     assert parentalGenosDict[record.POS]['REF'] == record.REF
        # except AssertionError:
        #     message = f'REF alleles for position {record.POS} don\'t match in founding SNPs and STITCH VCF'
        #     logging.warning(message)
        #     sys.exit(message)

        try:
            assert record.POS in parentalGenosDict.keys()
            logging.info(f'found {record.POS}')
        except AssertionError:
            message = f'position {record.POS} was not found in the founding SNPs'
            logging.warning(message)
            continue

        try:
            assert parentalGenosDict[record.POS]['REF'] == record.REF
        except AssertionError:
            message = f'REF alleles for position {record.POS} don\'t match in founding SNPs and STITCH VCF'
            logging.critical(message)
            sys.exit(message)

        try:
            assert parentalGenosDict[record.POS]['ALT'] == record.ALT[0]
        except AssertionError:
            message = f'ALT alleles for position {record.POS} don\'t match in founding SNPs and STITCH VCF'
            logging.critical(message)
            sys.exit(message)




        # parent1_geno = record.genotype(sample_names[0])['GT']
        # parent2_geno = record.genotype(sample_names[1])['GT']
        #
        # lenALT = len(record.ALT)
        #
        # if (parent1_geno != parent2_geno) and lenALT == 1 and record.num_unknown == 0:
        #     outputInfo = [record.CHROM, record.POS, record.REF, record.ALT[0], parent1_geno, parent2_geno]
        #     outputLine = '\t'.join(map(str, outputInfo))
        #     outputFile.write(outputLine+'\n')

if __name__ == "__main__":
    main()
