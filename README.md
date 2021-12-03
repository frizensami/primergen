# PrimerGen

## Goal
Generate as many primers as possible with these constraints:

1. Each primer has 20 nucleotides
1. Each primer must have a ratio of ((G + C) / (A + T + G + C)) between 0.45 and 0.55 ("balanced GC content")
1. Any two pairs of primers must have an edit distance >= 8

Metrics:

1. Total number of valid primers generated 
2. Primer generation speed (to allow regeneration if the constraints change)

## Previous work and algorithms
- 2009: **DeLOB** algorithm: Xu Q. , Schlabach M.R., Hannon G.J., Elledge S.J. (2009) Design of 240, 000 orthogonal 25mer DNA barcode probes. Proc. Natl. Acad. Sci. U S A, 106, 2289–2294. (https://www.pnas.org/content/pnas/106/7/2289.full.pdf)
  - DeLOB uses: Blast tool with the -F off option turned
- 2012: ***Schwartz et al**: “Accurate gene synthesis with tag-directed retrieval of sequence-verified DNA molecules” 
  - 2016: Bornholt et al: freeloads off Schwartz
- 2013: **TagGD** - super important, has actual GitHub code! Gene rates unique DNA sequences with edit distance and GC content (https://github.com/pelinakan/UBD - https://mybiosoftware.com/taggd-dna-tag-generation-demultiplexing.html) 
  - This seems pretty slow, for just 500 primers with min-8-edit distance, this took ages
- 2017: **Organick et al**, “Scaling up DNA data storage and random access retrieval”, very similar, using filtering and then BLAST 
- 2018: DeLOB modified for use in https://academic.oup.com/synbio/article/3/1/ysx008/4817474 (2018) to generate 20-mers 


## TODOs
- TODO: store primers found over time / iterations as a graph
- TODO: Re-run with further edit distance
  - With fixed time, no clear result difference. 10/9/9/10 primers random/gc/freq/noreroll. More PPI for last two but the advantage is not clear yet (4/7/9/8) * 10^-5.
- TODO: MOST IMPT Implement the network graph node-and-neighbour removal algorithm
  - MOST IMPT See if we can improve that by not just randomly removing nodes  (centrality etc)
- TODO: see if we can optimize for furthest possible edit distance
- TODO: Can we generate all balanced-GC strigs?
- TODO: make new strategy where sometimes we randomize completely
- TODO: store all the primer to primer distances in a file so we can plot a force directed graph of how far apart everything is


## DeLOB protocol (from https://www.pnas.org/content/pnas/106/7/2289.full.pdf)

Ten million 25mer oligo DNA sequences were generated
as candidates with the ‘‘makenucseq’’ program in the EMBOSS package (20).
These DNA sequences were sequentially fed into a restriction enzyme filter
which exclude sequences containing restrictive enzyme sites that are reserved
for library cloning (EcoR1, XhoI, BglII, MluI, AvrII, FseI, and MfeI), a Tm filter
based on the ‘‘nearest neighbor model’’ (21) to exclude sequences of Tm below
58 °C or above 68 °C, a GC composition filter to exclude sequences of GC below
40% or above 60%, and a repetitive sequence filter to exclude sequences
containing repetitive tracts (5 or longer single nucleotide repeats or 4 or
longer double nucleotides repeats). Candidates that passed all these filters
were compared to each other for sequence similarity using the BLAST program
with the ‘‘-F’’ option turned off. We defined 2 candidates to be orthogonal
to each other if they do not have stretches longer than 12 bases of HSPs
between them. On the basis of BLAST results, candidates were divided into 2
groups: those with no HSPs of 13 bases or longer to any other candidate
(orthogonal probes I), and those with longer than 12 bases HSPs to at least 1
of other candidates (nonorthogonal probes). For the latter group, we applied
a ‘‘network elimination’’ algorithm (see below) to obtain a subset of candidates that were orthogonal to each other (orthogonal probes II), and combine
with orthogonal probes I. These orthogonal probes were then fed into a
secondary structure filter, which was based on the ‘‘hybrid-ss’’ program in the
UNAFold package (12) to exclude probes that form intraprobe secondary
structures (self-folding energy  2 kJ/mol at 50 °C).

The Network Elimination Algorithm. We first constructed a network from all
nonorthogonal candidates. Each vertex in the network represented a candidate
and an edge represented the existence of a longer than 12-base HSP between the
2 connected candidates. We randomly chose 1 candidate and placed it in the
inclusion group (orthogonal probes II). Candidates that were connected to this
one were placed into the exclusion group. We then eliminated all candidates in
the exclusion group from the network, together with all edges incident to these
candidates. This selection-and-elimination procedure was then repeated on the
remaining network till all candidates were put into either of the 2 groups.
Candidates in the inclusion group were orthogonal to each other.

## Development
### Pre-requisites for building
1. `biopython` e.g., `conda install -c conda-forge biopython`
2. `editdistance` https://github.com/roy-ht/editdistance
3. `polyleven` https://github.com/fujimotos/polyleven (possibly faster levenshtein)
4. `networkx` - for graph algorithms


