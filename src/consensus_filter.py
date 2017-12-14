#!/usr/bin/env python
# -*- coding: utf-8 -*- \#
"""
@author = 'liangzb'
@date = '2017/12/14 0014'

"""

from Bio import SeqIO


def filter_few_reads(fastq_in, fastq_out, threshold=2):
    with open(fastq_in) as fp, open(fastq_out, 'w') as fp_out:
        for record in SeqIO.parse(fp, 'fastq'):
            umi, pos_reg, reads1, reads2 = record.id.split(':')
            if int(reads1) < threshold or int(reads2) < threshold:
                continue
            fp_out.write(record.format('fastq'))
