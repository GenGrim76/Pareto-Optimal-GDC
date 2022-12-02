# Pareto-Optimal-GDC: Pareto Optimal Genomic Dictionaries Compression

## Environment

The **Pareto-Optimal-GDC** java library runs over Linux/Ubuntu distros with java 11 pre-installed. It released under the terms of the GNU General Public License v3.0 (see *License.md*). It requires a *Java* compliant Virtual Machine (version >= 1.8).  It integrates libraries written in C++ and python, which is why the following pre-installed environments are required: 

- **python** version >= 3.7
- **cmake** version >= 3.25
- **gcc** version >= 11.0


## Library Dependencies

- **JavaFastPFOR** [1] library: A simple integer compression library in Java, released under the Apache License Version 2.0 http://www.apache.org/licenses/.

- **jUnit** Testing library: simple framework to write repeatable tests on the JVM, released under the Eclipse Public License 1.0 (EPL) https://www.eclipse.org/legal/epl-v10.html. 

- **jFreeChart** library: Java chart 100% open-source library that makes it easy for developers to display professional quality charts in their application, released under the terms of the GNU Lesser General Public Licence (LGPL) https://www.gnu.org/licenses/lgpl-3.0.html.

  #### (PLEASE NOTE: see the pom.xml file for the relative repositories)


## Additional Software used in this library

- aaa
- bbb
- ccc
- ddd
- eee




## Usage

The software is released as an executable **Pareto-Optimal-GDC.jar** java library, together with a whole set of integrated tools for the compression/decompression of k-mers, genomic sequences and integers, both generic and specialized, as further detailed in [2]. 

The Pareto-Optimal-GDC java library can be used to run experiments from the command line, using the following syntax:


`java -jar Pareto-Optimal-GDC.jar [input_parameters]`

where [input_parameters] as reported below:
1) datasets_directory
2) dataset_name 
3) experiment_directory 
4) debug [0/1] (1 for verbose mode)


## Quickstart
For a quick start, assuming both JVM, python and cmake are properly installed on your system, unzip the Pareto-Optimal-GDC 1.0 release, move on the directory `/bin`, and run the following command:

`java -jar Pareto-Optimal-GDC.jar ../data StaphAU ../experiments 1`

Execution of the previous command returns a set of files in '.txt' format. For the **Base Case Scenario**, and for each k-value used in this study, a special folder is created in which there are two files storing the performance in terms of compression and post-processing time obtained by the various compressors (standard and specialized) for the k-mers and relative frequencies.



## Datasets

**Pareto-Optimal-GDC.jar** library has been extensively tested by using the following genomic datasets:
- [Staphylococcus Aureus (3 MB)](https://www.ncbi.nlm.nih.gov/nuccore/NC_010079.1?report=fasta)
- [Human Chromosome 14 (109 MB)](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.14/)
- [Assembled Plants (4.8 GB)](http://afproject.org/media/genome/std/assembled/plants/dataset/assembled-plants.zip) (concatenation of the 14 plant genomes via **cat** linux command)


## References
[1]: Daniel Lemire and Leonid Boytsov. "Decoding billions of integers per second through vectorization." Software: Practice and Experience 45.1 (2015): 1-29.

[2]: paper reference
