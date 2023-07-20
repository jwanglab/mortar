Mortar builds a reference-based consensus sequence from raw read alignments.

It was designed specifically to work with nanopore sequencing of tiled multiplex PCR amplicons like those commonly used for SARS-CoV-2 and MPXV.
It also works equally well on shotgun genomic sequencing (for MPXV).
We haven't tried with Illumina sequence (either shotgun or amplicons), but there's no technical reason it wouldn't work well.

Overview
--------

1.  Alignments of reads exceeding the maximum are immediately excluded. This accounts for not-infrequent chimeric reads.
2.  If a BED file of multiplex PCR primers is provided, alignments are implicitly trimmed to include only the inter-primer region of the tile with the greatest overlap
3.  A simple count is taken of the allele (or deletion) and insertion at every reference locus
4.  At each site, the most commonly observed allele is assigned if it is not a deletion. Loci with total coverage below the designated threshold are reported as Ns.
5.  Deletions are handled according to the following heuristic rules:
    *  all sites with >50% deletion calls OR (deletion is the highest frequency call AND the allele matches the target locus [homopolymers]) are considered putative deletions; all other are assigned the highest-frequency nucleotide
    *  deletion blocks are composed of all contiguous deletions by the above rule
    *  a deletion block is accepted if the maximum deletion fraction among all contiguous deletion loci exceeds the given threshold (default: 0.7) OR have a total homopolymer count meeting the given threshold (default: 5)


Installation
------------

    git clone https://github.com/jwanglab/mortar
    cd mortar
    make
    ./mortar -h


Usage
-----

    Usage: mortar [options] <ref> <paf>
    Options:
      Reference FASTA file
      PAF alignment of reads to this reference, with --cs flag
      -b: BED file containing multiplex (tiled) PCR primers used
          If a BED file is included, primers will be identified and trimmed off before consensus
      -m: minimum depth below which reference sequence is masked with Ns [default: 20]
      -d: maximum depth to consider (additional will simply be ignored) [default: 1000]
      -r: maximum read length (longer reads [potentially chimeric] will simply be ignored) [default: 2700]
      -l: fraction of depth to call deletion (default: 0.7)
          the maximum fraction in a contiguous deletion must exceed this value to call a deletion
          all other sites within a deletion must just be the highest allele
      -c: output coverage per site to the given file
      -v: version
      -h: help

Example for MPXV using Chen primer set:

    minimap2 -cx map-ont --cs MT903340.fasta reads.fastq.gz > mpxv.paf
    ./src/mortar/mortar MT903340.fasta mpxv.paf -b Chen.bed -c covg.dat > mortar_consensus.fasta

To plot coverage:

    gnuplot -e 'set terminal png; set output "covg.png"; set autoscale; set yrange [*<0:2<*]; set logscale y 10; unset label; set xtic auto; set ytic auto; set title "Sequencing depth"; set xlabel "Position (nt)"; set ylabel "Depth"; set xtics rotate; plot "covg.dat" title "Depth" with lines; set terminal dumb; set output; replot'

BED file must be a tab-delimited file in this format:

    sequence_name	start	end	tile_name

Where the end positions are *inclusive* and the tile\_name must alternate matched primers that include the word "LEFT" and "RIGHT".

For example, the first two tiles for Chen's MPXV scheme:

    MT903340	356	378	MPXV_1_LEFT
    MT903340	2395	2417	MPXV_1_RIGHT
    MT903340	1439	1461	MPXV_2_LEFT
    MT903340	3353	3375	MPXV_2_RIGHT
    ...
