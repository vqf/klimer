# klimer
This set of applications illustrates the klimer algorithm for indexing of NGS reads. The applications are under development, but some of their capabilities are already available. The main tools are:

* buildMap - Accepts NGS reads in Fasta or Fastq formats and outputs an indexing file. The default output file has the name of the input file followed by the k11 extension. The length of the indexing k-mers is 11 by default, but can be changed with the -k option.
   Use: buildMap inputFile [output_file] [-k kmerLength]
* checkSeq - Looks for evidence for a sequence in an index file. Uses stdin for repeated querying. Use: checkSeq indexFile
* kmerHist - Builds a k-mer histogram for an arbitrary k value from an index file. Use: kmerHist indexFile [-k kmerLength]

# Installation and testing
