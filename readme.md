# a simple class for Unique Molecular Identifiers processing

requirements:

* Python3.6:  Python3 is OK, but I'd like to use [fstring](https://www.python.org/dev/peps/pep-0498/) in Python3.6, if you're using Python3 but under 3.6, feel free to replace the `fstring` statement with `format` statement.
* path.py
* biopython: BioPython is for fastq files processing.
* pysam: standard library for bam file processing.
* fastinterval: genomic interval processing.

# UMI pre-processing

it only works with **sorted** bam file, with UMI at the end of reads id:

```text
M03074:71:000000000-B49DJ:1:2103:19133:13331:TTCGCACGA
M03074:71:000000000-B49DJ:1:2105:4089:17175:TTCGCACGA
```

for UMI pre-processing, here are some choices:

## 1. use `bcl2fastq`

if your data is from illumina platform, just use `bcl2fastq`.

Before you doing this, you must make sure the UMI cycles is in either read1 or read2 cycles, you can change **Reads** section in **RunInfo.xml**:

for example, if your UMI is in INDEX2, with 10 bases, your RunInfo.xml will looks like:

```text
<Read Number="1" NumCycles="150" IsIndexedRead="N" />
<Read Number="2" NumCycles="8" IsIndexedRead="Y" />
<Read Number="3" NumCycles="10" IsIndexedRead="Y" />
<Read Number="4" NumCycles="150" IsIndexedRead="N" />
```

since bcl2fastq can not handle the UMI in the index region, you must combine the UMI cycles to reads2 cycles, like:

```text
<Read Number="1" NumCycles="150" IsIndexedRead="N" />
<Read Number="2" NumCycles="8" IsIndexedRead="Y" />
<Read Number="4" NumCycles="160" IsIndexedRead="N" />
```

now, bcl2fastq will take the UMI cycles for Reads2.

set settings section in your *SampleSheet.csv*:

```text
[Settings]
Read1StartFromCycle,1
Read1EndWithCycle,150
Read2StartFromCycle,1
Read2EndWithCycle,160
Read2UMIStartFromCycle,1
Read2UMILength,10
TrimUMI,1
```

make sure your csv file with same number of columns.

## 2. use `fastp`

[fastp](https://github.com/OpenGene/fastp) is a fast all-in-one tool for preprocessing FastQ files.

It can extract UMI too.

# Pre-Mapping and sorting

using `bwa` or other aligner, produce a sorted bam file, for example:

```bash
bwa mem -t 8 -M hg19.fasta read1.fastq read2.fastq | \
samtools view -@ 8 -S -b | \
samtools sort -@ 8 -o sample_name.sorted.bam 
```

# make consensus

```python
from consensus_maker import ConsensusWorker

worker = ConsensusWorker('sample_name.sorted.bam', 'consensus.1.fastq', 'consensus.2.fastq', 
                         bed_file='target.bed')
worker.output_pe_reads()
```

# TODO list

**TODO**:

1. make a simple script;
2. UMI analysis and report;
3. finish testing;
4. support for `a-b` family and `b-a` family.

welcome to pull a request. 
