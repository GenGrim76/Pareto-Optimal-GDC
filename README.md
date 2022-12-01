# Compressed_kmer_Dictionaries
Disk Storage of Compressed k-mer Dictionaries, with or without Random Access in Main Memory.


## Environment
-
-
-


## Libraries Dependencies

- JavaFastPFOR library: A simple integer compression library in Java, released under the Apache License Version 2.0 http://www.apache.org/licenses/.
- junit Testing Library: simple framework to write repeatable tests on the JVM. 
- jcommander Command Line Interface Library: Command line parsing framework for Java.
- jFreeChart Library: Java chart library that makes it easy for developers to display professional quality charts in their application.

  #### (PLEASE NOTE: see the pom.xml file for the different library repositories)




## Datasets

TradeOff jar library has been extensively tested by using the following genomic datasets:
- [Staphylococcus Aureus (3 MB)](https://www.ncbi.nlm.nih.gov/nuccore/NC_010079.1?report=fasta)
- [Human Chromosome 14 (109 MB)](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.14/)
- [Assembled Plants (4.8 GB)](http://afproject.org/media/genome/std/assembled/plants/dataset/assembled-plants.zip) (concatenation of the 14 plant genomes via cat linux command)
