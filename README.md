# klimer
This set of applications illustrates the klimer algorithm for indexing of NGS reads. The applications are under development, but some of their capabilities are already available. The main tools are:

* buildMap - Accepts NGS reads in Fasta or Fastq formats and outputs an indexing file. The default output file has the name of the input file followed by the k11 extension. The length of the indexing k-mers is 11 by default, but can be changed with the -k option.
   Use: buildMap inputFile [output_file] [-k kmerLength]
* checkSeq - Looks for evidence for a sequence in an index file. Uses stdin for repeated querying. Use: checkSeq indexFile
* kmerHist - Builds a k-mer histogram for an arbitrary k value from an index file. Use: kmerHist indexFile [-k kmerLength]

# Installation and testing
The folder contains the project files from Code:Blocks. In Linux, the script compile.sh followed by the name of a c file should compile it. The folders Win32 and LinuxBin contain pre-compiled binaries for both systems. Windows binaries cannot read compressed files, but Linux binaries can read gz files. The folder contains a few E. coli reads for quick tests. For more thorough testing, you can download a medium-size Fastq file, like ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR306/002/SRR3068732/SRR3068732_1.fastq.gz

>./builMap SRR3068732_1.fastq.gz

This may take up to one hour depending on the system. The peak RAM use is around 9 Gb. The result file will be called SRR3068732_1.fastq.gz.k11

>./searchSeq SRR3068732_1.fastq.gz.k11

The program will prompt for sequences. Since SRR3068732_1.fastq.gz contains E. coli sequences, here are some examples, taken from the genomic sequence, that should work:

*AGTGGGCCGGAGTTCACGGTTGCGATACTCGGTGAAGAAATTTTACCGTCAGTACG

*GAAAAACAGTTTTAGATAATAAGGAATATCTCAATTATTGAACATTTAGTGCGAATTATTTAGTACAAAAAAGCGGCGTTAGGTGATCTTTCCCTGGCCTC

*AATTACCGGTCAAGGGCATTTCCAATCTGAATAATATGGCAATGTTCAGCGTTTCCGGCCCGGGGATGAAAGGAATGGTCGGCATGGCGGCGCGCGTCTTTGCTGCAATGTCACGCGCCCGTATTTCCGTGGTGCTGATTACGCAATCATCTTCCGAATACAGTATCAGTTTCTGCGTTCCGCAAAGCGACTGTGTGCGAGCTGAACGGGCAATGCAGGAAGAGTTCTACCTGGAACTGAAAGAAGGCTTACTGGAGCCGCTGGCGGTGACGGAACGGCTGGCCATTATCTCGGTGGTAGGTGATGGTATGCGCACCTTGCGTGGGATCTCGGCGAAATTCTTTGCCGCGCTGGCCCGCGCCAATATCAACATTGTC

*TGCTGTAACTGCTGCTGATATGTGTTGAGCGATCGGCAAACGTACCAACCCACCCAACAGAGAAGCCACTGTGAGCGATCGGAATATTCAGAGTGCTGGTAACAGTATCCGGGTTAATGCTGGAGATGTATTCGCCGGTATCGGTGTCTTTGCCGCGGGTACGGTTATAGGCCACATCAAGGCTAAACAGATCAGTGGTATATTTCGTCATCACATCCCAGCCCCAGATTTTGGCGTTCGGGACGTTATACGACATAGTCGTCGCCGCCGCGAAATCGACGGTCGTGGAGATGTAATCCTTCGCTTTGGTATCAAAGTAGCTGGCTTTAAATTCCAGAGCATCATTGGACAACATCAGGTCATCAAAACGCAGCCCAAAACCGTACTCCTGAGTTTCGTTAGTTTCCGGACGTAAGTTCGGGTTTGGCACCCAATAGTTGGTATAGAAGCGACCAATCGAGAAGTGCTTAGAATCGTTATACATTTCGCCCATCGTCGGGGCGCGGAATGCCTGGGCATATGAGCCAAATAACATCAGCCAGTTAGTCGGATTGATAGTCATCCCCGCACGAGATGACCATTTGTCGGCATCAACATCTTTGTAACCGTCACTGCTACCGCGATAACTGTCATAGCGGGTTCCGCCAAGCAGGGTAATCGGCAGATCGCGTAAGGTGATCTCATCCTGTAGCCAGCCGGAGCTAAAATCGATTTTTGCTTGCGGGAAGCCCGTCGTCGCGCCGCCCGGATGTTGTTCCTGACGATAATACTCACCGCCATATGTCAGTAAGTGAGAAGCGAAACTGTCGGCAAAGAGAGTGGAACGGTTCTCCAGCCTGGCTCCTTTTGTTATCTGTTCACGATACTCGCCGGAACTCCCCGTGTTTTGCGCATTAATACGGACTTCCGACCAATAAATTTTTGCATCTGCATTTAACCAGTCGTTGCCCTGCGGGGCGAGTTTATAAGAAAGCTGCGCATCGCGTTGAATTGTTGAACGATCAACCATCGG


>./kmerHist SRR3068732_1.fastq.gz.k11 -k 25

Will print a histogram of 25-mers in the sample. The value for k is arbitrary, and can be larger than the input sequences. For values higher than 100, this may take several minutes.