# Compressed_kmer_Dictionaries
Disk Storage of Compressed k-mer Dictionaries, with or without Random Access in Main Memory.


## Environment
TradeOff jar library runs over Linux/Ubuntu distros with java 11 pre-installed. It requires a *Java* compliant Virtual Machine (>= 1.8).  It integrates libraries written in C++ and python, which is why the following pre-installed environments are required: 

- python version >= 3.7
- cmake version >= 3.25
- gcc version >= 11.0


## Libraries Dependencies

- JavaFastPFOR library: A simple integer compression library in Java, released under the Apache License Version 2.0 http://www.apache.org/licenses/.

- junit Testing Library: simple framework to write repeatable tests on the JVM, released under the Eclipse Public License 1.0 (EPL) https://www.eclipse.org/legal/epl-v10.html. 

- jFreeChart Library: Java chart 100% open-source library that makes it easy for developers to display professional quality charts in their application, released under the terms of the GNU Lesser General Public Licence (LGPL) https://www.gnu.org/licenses/lgpl-3.0.html.

  #### (PLEASE NOTE: see the pom.xml file for the different library repositories)


## Usage
The software is released as a single executable jar file, **TradeOff.jar**, that can be used to run experiments from the command line, using the following syntax:

`java -jar TradeOff.jar [input_parameters]`

where [input_parameters] as reported below:
1) datasets_directory
2) dataset_name 
3) experiment_directory 
4) DEBUG [0/1] (verbose mode)


## Datasets

**TradeOff.jar** library has been extensively tested by using the following genomic datasets:
- [Staphylococcus Aureus (3 MB)](https://www.ncbi.nlm.nih.gov/nuccore/NC_010079.1?report=fasta)
- [Human Chromosome 14 (109 MB)](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.14/)
- [Assembled Plants (4.8 GB)](http://afproject.org/media/genome/std/assembled/plants/dataset/assembled-plants.zip) (concatenation of the 14 plant genomes via cat linux command)
