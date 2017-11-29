import sys

from UMI.src.consensus_maker import ConsensusWorker

input_bam_file, read1_file, read2_file, stat_json = sys.argv[1:]
worker = ConsensusWorker(input_bam_file)
worker.output_pe_reads(read1_file, read2_file)
worker.output_stats(stat_json)
