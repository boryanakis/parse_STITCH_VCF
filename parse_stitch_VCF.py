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
        logging.debug(elements)
        if elements[0] == chrom:
            tempDict = OrderedDict(zip(header, elements))
            logging.debug(f'{tempDict}')

            genoDict[tempDict['POS']] = {k:tempDict[k] for k in header[2:]}
    logging.debug(f'\n{pprint.pformat(genoDict)}')

    try:
        assert len(genoDict.items()) > 0
    except AssertionError:
        message = f'There are no entries for chromosome {chrom} in the parental genotypes file.'
        logging.critical(message)
        sys.exit(message)

    return genoDict

def verifyOutputGenoIntegrity(genosDict, samples):

    nPos = len(genosDict['positions'])
    nBad = 0
    for sample in samples:
        try:
            assert len(genosDict[sample]) == nPos
        except AssertionError:
            nBad += 1
            message = f'sample [{sample}] has [{len(genosDict[sample])}] genos. Expected: [{nPos}]'
            logging.warning(message)

    if nBad == 0:
        result = True
    else:
        result = False

    return result


def main():
    """

    :return:
    """

    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-i', '--inVCF',
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
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='print more information')
    args = parser.parse_args()

    logger = logging.getLogger('root')
    FORMAT = "[%(filename)s:%(lineno)4s - %(funcName)20s() ] %(levelname)10s - %(message)s"
    logging.basicConfig(level=logging.WARNING, format=FORMAT)

    if args.verbose:
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

    logging.info(f"parental genotypes have been stored in memory")

    genotypeTranslationDictProper = {'0/0': '0', '0/1': '1', '1/1': '2', './.':'NA'}
    genotypeTranslationDictInv = {'0/0': '2', '0/1': '1', '1/1': '0', './.':'NA'}

    logging.info(f"opening the input VCF file...")
    vcfFileInput = openIOFile(args.inVCF)
    vcf_reader = vcf.Reader(vcfFileInput)

    nSamples = len(vcf_reader.samples)
    logging.info(f"calculated number of samples: [{nSamples}]")

    # prepare the data structures for the output file
    outputHeader = ['sample']

    outputGenosDict = OrderedDict()
    outputGenosDict['positions'] = []
    for sample in vcf_reader.samples:
        outputGenosDict[sample] = []

    logging.info(f"iterating through the VCF records")

    i = 0
    for record in vcf_reader:
        i += 1
        if i % 10000 == 0:
            logging.info(f"processing record # {i}")
        # print(record.CHROM, record.POS, record.num_called, record.num_unknown)

        # ensure that the current VCF position is in the parental genotypes table
        try:
            assert str(record.POS) in parentalGenosDict.keys()
            logging.debug(f'found {record.POS}')
        except AssertionError:
            message = f'position {record.POS} was not found in the founding SNPs'
            logging.warning(message)
            continue

        # check whether the REF allele matches between the VCF and the parental genos table
        try:
            assert parentalGenosDict[str(record.POS)]['REF'] == record.REF
        except AssertionError:
            message = f'REF alleles for position {record.POS} don\'t match in founding SNPs and STITCH VCF'
            logging.critical(message)
            sys.exit(message)

        # check whether the ALT allele matches between the VCF and the parental genos table
        try:
            assert parentalGenosDict[str(record.POS)]['ALT'] == record.ALT[0]
        except AssertionError:
            message = f'ALT alleles for position {record.POS} don\'t match in founding SNPs and STITCH VCF'
            logging.critical(message)
            sys.exit(message)

        # determine which translation dictionary will be used
        if parentalGenosDict[str(record.POS)][parents[0]] == '0/0' and parentalGenosDict[str(record.POS)][parents[1]] == '1/1':
            translateGeno = genotypeTranslationDictProper.copy()
        elif parentalGenosDict[str(record.POS)][parents[0]] == '1/1' and parentalGenosDict[str(record.POS)][parents[1]] == '0/0':
            translateGeno = genotypeTranslationDictInv.copy()
        else:
            message = f'An unexpected genotype combination was encountered in the parents at position [{record.POS}]'
            logging.warn(message)
            warnings.warn(message, Warning)
            continue

        outputGenosDict['positions'].append(str(record.POS))
        # iterate through the samples in the VCF
        # re-code genos as 0, 1, or 2
        for sample in vcf_reader.samples:
            trGeno = translateGeno[record.genotype(sample)['GT']]
            outputGenosDict[sample].append(trGeno)

    logging.debug(pprint.pformat(outputGenosDict))


    try:
        assert verifyOutputGenoIntegrity(outputGenosDict, vcf_reader.samples)
    except AssertionError:
        message = f'the output genotype dictionary is not correct'
        logging.critical(message)
        sys.exit(message)

    outputFN = f'{args.prefix}.{args.chrom}.genos.tsv'
    outputFile = openIOFile(outputFN, args.outputDir, 'w')
    outputHeader += outputGenosDict['positions']
    outputFile.write(','.join(map(str, outputHeader)) + '\n')

    for sample in vcf_reader.samples:
        outputLine = [sample] + outputGenosDict[sample]
        outputFile.write(','.join(map(str, outputLine)) + '\n')







if __name__ == "__main__":
    main()
