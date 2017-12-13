import argparse

from path import Path
from src.consensus_analysis import plot_reads_distribution
from src.consensus_maker import ConsensusWorker


def parse_args():
    """
    dest must be the same with the argument of ConsensusWorker and ConsensusMaker
    :return:
    """
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-i', '--bam_in', dest='bam_in',
                        metavar='FILE', type=Path, required=True,
                        help="input bam file")
    parser.add_argument('-q1', '--reads1', dest='reads1',
                        metavar='FILE', type=Path, required=True,
                        help="out reads1.fastq file")
    parser.add_argument('-q2', '--reads2', dest='reads2',
                        metavar='FILE', type=Path, required=True,
                        help="out reads2.fastq file")
    parser.add_argument('-s', '--stat_file', dest='stat_file',
                        metavar='FILE', type=Path, required=True,
                        help="out stat json file")
    parser.add_argument('-b', '--bed_file', dest='bed_file',
                        metavar='FILE', type=Path, default=None,
                        help="target intervals")
    parser.add_argument('-t', '--threads', dest='threads',
                        metavar='INT', type=int, default=8,
                        help="set threads number, [8]")
    parser.add_argument('--report', dest='report',
                        metavar='FILE', type=Path, default=None,
                        help="figure report file")
    parser.add_argument('--flank_size', dest='flank_size', metavar='INT', type=int, default=20,
                        help="flank size, [20]")
    parser.add_argument('--min_reads', dest='min_reads', metavar='INT', type=int, default=1,
                        help="minimum count of reads to make a consensus, [1] ")
    parser.add_argument('--min_qual', dest='min_qual', metavar='INT', type=int, default=10,
                        help="minimum phred quality of base can be used, [10]")
    parser.add_argument('--qual_n', dest='qual_n', metavar='INT', type=int, default=2,
                        help="N bases quality if make consensus base failed, [2]")
    parser.add_argument('--max_n_count', dest='max_n_count', metavar='INT', type=int, default=3,
                        help="maximum N bases in consensus reads, [3]")
    parser.add_argument('--max_n_ratio', dest='max_n_ratio', metavar='FLOAT', type=float, default=0.03,
                        help="maximum N ratio in consensus reads, [0.03]")
    parser.add_argument('--min_freq', dest='min_freq', metavar='FLOAT', type=float, default=0.6,
                        help="min frequency to make a consensus, [0.6]")
    args = parser.parse_args()
    assert args.min_freq >= 0.5, "--min_freq must be higher than 0.5"
    return args


args = parse_args()
worker = ConsensusWorker(**vars(args))
# worker.async_output_pe_reads(args.reads1, args.reads2)
worker.output_pe_reads(args.reads1, args.reads2)
worker.output_stats(args.stat_file)
if args.report:
    plot_reads_distribution(args.reads1, args.report)
