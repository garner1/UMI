#!/usr/bin/env python
# -*- coding: utf-8 -*- \#
"""
@author = 'liangzb'
@date = '2017/11/29 0029'

analysis results from consensus_maker
"""

from Bio import SeqIO
from matplotlib import pyplot as plt

plt.style.use('seaborn-paper')


def get_reads_count(fastq_file):
    reads1 = []
    reads2 = []
    with open(fastq_file) as fp:
        for record in SeqIO.parse(fp, 'fastq'):
            umi, r1, r2 = record.id.split(':')
            reads1.append(int(r1))
            reads2.append(int(r2))
    return reads1, reads2


def plot_reads_distribution(fastq_file, out_file, x_max=8):
    def subplot_distribution(ax, array, title):
        ax.set_title(title)
        ax.set_xlabel('#reads used')
        ax.set_ylabel('#consensus read')
        ax.set_xlim([1, x_max])
        hist = ax.hist(array)
        rects = hist[2]
        value = hist[0]
        sum_value = sum(value)
        for ind, rect in enumerate(rects):
            ax.text(rect.get_x() + rect.get_width() / 4.,
                    rect.get_height(),
                    f'{value[ind]/sum_value*100:2.2f}%')

    fig = plt.figure(figsize=(20, 9))
    fig.subplots_adjust(wspace=0.25, hspace=0.20, top=0.85, bottom=0.05)
    fig.suptitle("reads used when making a consensus")
    reads1, reads2 = get_reads_count(fastq_file)

    ax1 = fig.add_subplot(1, 2, 1)
    subplot_distribution(ax1, reads1, 'Reads1')
    ax2 = fig.add_subplot(1, 2, 2)
    subplot_distribution(ax2, reads2, 'Reads2')

    fig.savefig(out_file)
