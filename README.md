Library Preparation
`—TODO—`

Sequencing and Identification of Somatic Mutations
===

DCS Tag Filtering and Concatentation to Header
---
The library of amplified fragments from the harvested livers of four mice exposed to AFB1 where sequenced using the Illumina NextSeq platform as paired-end 150 bp reads. Reads were then processed through `tag_to_header.py` in order to extract and concatenate the 8 nt duplex consensus sequence (DCS) tag. Reads were filtered out if their 5 nt Illumina spacer was missing. Reads were also filtered out if their DCS tag contained any non-alphanumeric characters, any Ns, or any homopolymeric runs of at least 8 nt long. `tag_to_header.py` can also be used to filter out any read without a known spacer if that string is supplied with the `--spacer <string>` command.

```bash
python ${DSpath}/tag_to_header.py \
	--infile1 $read1in \
	--infile2 $read2in \
	--outfile1 ${runIdentifier}.seq1.fq.smi \
	--outfile2 ${runIdentifier}.seq2.fq.smi \
	--blength 8 \
	--slength 5 \
	--poly 8 \
	--read_out 1000000
```

Alignment to Genome and Sorting
---
A modified EG10 genome (48417 nt) was constructed and indexed with the elimination of two germline mutations in the lambda transgene. Paired-end sequences were then aligned individually to this genome using the burrows-wheeler aligner (v0.6.1).

```bash
bwa aln $alignRef ${runIdentifier}.seq1.fq.smi > ${runIdentifier}.seq1.aln
bwa aln $alignRef ${runIdentifier}.seq2.fq.smi > ${runIdentifier}.seq2.aln
```

The two sets of paired-end sequences were then transformed using the burrows-wheeler command sampe to find the pairs of alignments which select the best hits according to orientations, distances between reads, etc. (2 smi files > 1 sam file) with the Smith-Waterman algorithm disabled for the unmapped mate.

```bash
bwa sampe \
	-s $alignRef \
	${runIdentifier}.seq1.aln \
	${runIdentifier}.seq2.aln \
	${runIdentifier}.seq1.fq.smi \
	${runIdentifier}.seq2.fq.smi \
	> ${runIdentifier}.pe.sam
```

Collapse SSCS and Analyze Tagstats
---
The paired-end aligned reads were then sorted with Samtools (v0.1.19) *view* as sam files `-S` and outputted in an uncompressed `-u` bam `-b` file.

```bash
samtools view -Sbu ${runIdentifier}.pe.sam \
	| samtools sort - ${runIdentifier}.pe.sort
```

Single strand consensus sequences (SSCS) were then collapsed using the Loeb Lab’s script `ConsensusMaker.py` **<REMAKE IN HOUSE>** using the `os` filter set and `dpm` read types settings.

```bash
python ${DSpath}/ConsensusMaker.py \
	--infile ${runIdentifier}.pe.sort.bam \
	--tagfile ${runIdentifier}.pe.tagcounts \
	--outfile ${runIdentifier}.sscs.bam \
	--minmem 3 \
	--maxmem 1000 \
	--readlength 137 \
	--cutoff 0.7 \
	--Ncutoff 0.3 \
	--read_type 'dpm' \
	--filt 'os' \
	--isize -1
```

Tagstats were examined to ensure an optimal family size distribution as shown **<GENERATE FIGURES FOR EACH MOUSE>**. Loeb recommends a family size of no more than 16 for efficient use of information. SSCS were finally sorted and indexed using Samtools view (-bu) and index commands.

```bash
samtools view -bu ${runIdentifier}.sscs.bam \
	| samtools sort - ${runIdentifier}.sscs.sort.bam

samtools index ${runIdentifier}.sscs.sort.bam
```

Create Duplex Consensus Sequences and Align back to Reference
---
DCSs were created using Loeb’s script `DuplexMaker.py` **<REMAKE IN HOUSE>**. The `DuplexMaker.py` algorithm searches for complimentary α/β tags and constructs a matching sequence only where bases match exactly. Any position with a mismatch is assigned an *N*.

```bash
python ${DSpath}/DuplexMaker.py \
	--infile ${runIdentifier}.sscs.sort.bam \
	--outfile ${runIdentifier}.dcs.bam \
	--Ncutoff 0.5 \
	--readlength 137 \
	--barcode_length 8
```

DCSs are then aligned back the reference genome, sorted, and indexed as described previously.

```bash
bwa aln $alignRef ${runIdentifier}.dcs.r1.fq > ${runIdentifier}.dcs.r1.aln
bwa aln $alignRef ${runIdentifier}.dcs.r2.fq > ${runIdentifier}.dcs.r2.aln

bwa sampe \
	-s $alignRef ${runIdentifier}.dcs.r1.aln \
	${runIdentifier}.dcs.r2.aln ${runIdentifier}.dcs.r1.fq \
	${runIdentifier}.dcs.r2.fq \
	> ${runIdentifier}.dcs.sam

samtools view -Sbu ${runIdentifier}.dcs.sam \
	| samtools sort - ${runIdentifier}.dcs.aln.sort

samtools index ${runIdentifier}.dcs.aln.sort.bam
```

Post-processing of Reads and Pileup Generation
---
DCSs are filtered using Samtools *view* `-F 4` which does not output any unmapped alignments. The output file is a compressed BAM file `-b`.

```bash
samtools view \
	-F 4 \
	-b ${runIdentifier}.dcs.aln.sort.bam \
	> ${runIdentifier}.dcs.filt.bam
```

All reads are then grouped using Picard’s (v1.119) *AddOrReplaceReadGroups.jar* command and indexed using Samtools *index*.

```bash
java -jar $picard/AddOrReplaceReadGroups.jar \
	INPUT=${runIdentifier}.dcs.filt.bam \
	OUTPUT=${runIdentifier}.dcs.filt.readgroups.bam \
	RGLB=UW \
	RGPL=Illumina \
	RGPU=ATATAT \
	RGSM=default

samtools index ${runIdentifier}.dcs.filt.readgroups.bam
```

Reads are then realigned around indels using the Genome Analysis Toolkit (GATK) (v3.30) *RealignerTargetCreater* and *IndelRealigner* functions.

```bash
java -Xmx2g -jar $gaTK/GenomeAnalysisTK.jar \
	-T RealignerTargetCreator \
	-R $alignRef \
	-I ${runIdentifier}.dcs.filt.readgroups.bam \
	-o ${runIdentifier}.dcs.filt.readgroups.intervals

java -Xmx2g -jar $gaTK/GenomeAnalysisTK.jar \
	-T IndelRealigner \
	-R $alignRef \
	-I ${runIdentifier}.dcs.filt.readgroups.bam \
	-targetIntervals ${runIdentifier}.dcs.filt.readgroups.intervals \
	-o ${runIdentifier}.dcs.filt.readgroups.realign.bam
```
Reads are then clipped using the GATK *ClipReads* function based on an analysis of base substitutions, Ns, and indels at each position in the read **<FIGURES FOR EACH>**. For this analysis we use settings which clip a total of 24 nt `--cyclesToTrim “1-4,117-137` `--clipRepresentation SOFTCLIP_BASES`. This reduces our reads to an effective 112 nt length.

```bash
java -Xmx2g -jar $gaTK/GenomeAnalysisTK.jar \
	-T ClipReads \
	-I ${runIdentifier}.dcs.filt.readgroups.realign.bam \
	-o ${runIdentifier}.dcs.filt.readgroups.clipped.bam \
	-R $alignRef \
	--cyclesToTrim "1-4,117-137" \
	--clipRepresentation SOFTCLIP_BASES
```

We then processed the final BAM file with the Samtools command *mpileup* with per-base alignment quality disabled `-B`, anomalous read pairs counted `-A` and a max depth of 500000 `-d 500000`.

```bash
samtools mpileup \
	-B \
	-A \
	-d 500000 \
	-f $alignRef \
	${runIdentifier}.dcs.filt.readgroups.clipped.bam \
	> ${runIdentifier}.dcs.pileup
```

Analysis of Pileup File
---
We designed the Python 3 library Siglib to handle the analysis of pileup files for spectra generation and comparison. First post-processing is achieved by parsing the pileup into a `siglib.pileup` object which allows for a descriptive view of pileup and read specific statistics. If `unique=True` then all validated somatic mutations at a single position will be listed only once in the resulting mutpos file.

> Note: Listing `unique=True` does not affect pileup and read specific statistics.

```python
from siglib.pileup import pileup

mpileup = pileup('${runIdentifier}.dcs.pileup', id=${runIdentifier})
mpileup.parse('${runIdentifier}.dcs.d100-c0.0-C0.01.mutpos',
              min_depth=100,
              clonality=(0.0, 0.01),
              unique=False)
mpileup.plot_statistics()
mpileup.report()
```
```
OUT [1]:

ID: ${runIdentifier}
Name:
Description:
Minimum Depth:    100
Clonality Filter: 0 - 0.01

                      COUNT    FREQUENCY   LOWER_CONF   UPPER_CONF
A's Sequenced:     16911023
       A to C:            0   0.0000e+00   0.0000e+00   2.2716e-07
       A to G:            0   0.0000e+00   0.0000e+00   2.2716e-07
       A to T:            3   1.7740e-07   6.0332e-08   5.2162e-07
C's Sequenced:     17588728
       C to A:            4   2.2742e-07   8.8439e-08   5.8480e-07
       C to G:            1   5.6855e-08   1.0036e-08   3.2208e-07
       C to T:           16   9.0967e-07   5.5996e-07   1.4778e-06
G's Sequenced:     18086668
       G to A:            8   4.4231e-07   2.2413e-07   8.7289e-07
       G to C:            0   0.0000e+00   0.0000e+00   2.1239e-07
       G to T:            4   2.2116e-07   8.6004e-08   5.6870e-07
T's Sequenced:     16248703
       T to A:            3   1.8463e-07   6.2791e-08   5.4289e-07
       T to C:            4   2.4617e-07   9.5732e-08   6.3303e-07
       T to G:            0   0.0000e+00   0.0000e+00   2.3642e-07

Nt Sequenced:      68835122
Pt Substitutions:        43   6.2468e-07   4.6380e-07   8.4137e-07

+1 insertions:            8   1.1622e-07   5.8891e-08   2.2935e-07
-1 deletions:             1   1.4527e-08   2.5645e-09   8.2297e-08

Total Insertions:         8   1.1622e-07   5.8891e-08   2.2935e-07
Total Deletions:          1   1.4527e-08   2.5645e-09   8.2297e-08
```
