# Pareto-Optimal-GDC: Pareto Optimal Genomic Dictionaries Compression

## Environment

The **Pareto-Optimal-GDC** java library runs over Linux/Ubuntu distros with java 11 pre-installed. It released under the terms of the GNU General Public License v3.0 (see *LICENSE.md*). It requires a *Java* compliant Virtual Machine (1.8+).  It integrates libraries written in C++ and python, which is why the following pre-installed environments are required: 

- **python** 3.7+
- **CMake** 3.25+
- **gcc** 11.0+


## Library Dependencies

- **JavaFastPFOR** [1] library: A simple integer compression library in Java, released under the Apache License Version 2.0 http://www.apache.org/licenses/.

- **jUnit** Testing library: simple framework to write repeatable tests on the JVM, released under the Eclipse Public License 1.0 (EPL) https://www.eclipse.org/legal/epl-v10.html. 

- **jFreeChart** library: Java chart 100% open-source library that makes it easy for developers to display professional quality charts in their application, released under the terms of the GNU Lesser General Public Licence (LGPL) https://www.gnu.org/licenses/lgpl-3.0.html.

  ### (PLEASE NOTE: see the pom.xml file for the relative repositories)


## Additional Software used in this library

- **DSK**: k-mer counter for reads or genomes. [https://gatb.inria.fr/software/dsk/] [2]
- **KMC3** disk-based program for counting k-mers from (possibly gzipped) FASTQ/FASTA files. [https://github.com/refresh-bio/KMC] [3]
- **ESSCompress** A tool to compress a set of k-mers represented in FASTA/FASTQ/KFF file(s). [https://github.com/medvedevgroup/ESSCompress] [4]
- **FM-index** [5] (as implemented in the **SDSL-Lite** library [6])
- **BCSF** (implemented in the **locom** library) [7]
- **bzip2** [8]
- **lz4** [9]
- **Zstandard** (**zstd**, for short) [10]
- **MFCompress** (**MFC**, for short) [11]
- **SPRING** [12]
- **BIC** (as implemented in the **JavaFastPFOR** library)
- **Opt-PFOR** (as implemented in the **JavaFastPFOR** library)



## Usage

The software is released as an executable **Pareto-Optimal-GDC.jar** java library, together with a whole set of integrated tools for the compression/decompression of k-mers, genomic and integers sequences, both generic and specialized, as further detailed in [XXXXXXXXXXXXXXXXX]. 

The Pareto-Optimal-GDC java library can be used to run experiments from the command line, using the following syntax:


`java -jar Pareto-Optimal-GDC.jar [input_parameters]`

where [input_parameters] as reported below:
1) datasets_directory
2) dataset_name 
3) experiment_directory 
4) debug [0/1] (1 for verbose mode)


## Quickstart

**PLEASE NOTE:** Depending on the size of the input genomic dataset, in addition to consuming a lot of time, use of the library requires a lot of disk space for the creation of temporary files, even on the order of several TB, as in the case of *Assembled Plants* Genome.

For a quick start, assuming both **JVM**, **python** and **CMake** are properly installed on your system, unzip the Pareto-Optimal-GDC 1.0 release, move on the directory `/bin`, and run the following command:

`java -jar Pareto-Optimal-GDC.jar ../data StaphAU ../experiments 1`


Execution of the above command returns a series of files in '.txt' format. For the **Base Case Scenario**, a special folder is created in which there are two folders: **CD-NRAM** and **SD-RAM**. In each of these there are two text files that respectively store the performance in terms of compression and post-processing time obtained by the various compressors (standard and specialised) on the k-mer dictionaries and on their relative frequencies, respectively for the "**Compressed on Disk - No Random Access in Main Memory**" and "**Succinct on Disk - Random Access in Main Memory**" Scenarios. For each value of k used in this study, a folder is created in which the results in terms of compression and post-processing time are similarly reported for each of the **DP0**, **DP1**, **DP2** and **DP3** Cases.



## Datasets

**Pareto-Optimal-GDC.jar** library has been extensively tested by using the following genomic datasets:
- [*Staphylococcus Aureus* (3 MB)](https://www.ncbi.nlm.nih.gov/nuccore/NC_010079.1?report=fasta)
- [*Human Chromosome 14* (109 MB)](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.14/)
- [*Assembled Plants* (4.8 GB)](http://afproject.org/media/genome/std/assembled/plants/dataset/assembled-plants.zip) (concatenation of the 14 plant genomes via **cat** linux command)


## References
[1]: Daniel Lemire and Leonid Boytsov. "Decoding billions of integers per second through vectorization." Software: Practice and Experience 45.1 (2015): 1-29.

[2]: Guillaume Rizk, Dominique Lavenier and Rayan Chikhi. DSK: k-mer counting with very low memory usage, Bioinformatics, 29(5):652-653 (2013). doi:10.1093/bioinformatics/btt020 .

[3]:  Marek Kokot, Maciej Długosz, Sebastian Deorowicz, KMC 3: counting and manipulating k-mer statistics, Bioinformatics, Volume 33, Issue 17, 01 September 2017, Pages 2759–2761, https://doi.org/10.1093/bioinformatics/btx304 .

[4]: Amatur Rahman, Rayan Chikhi and Paul Medvedev, Disk compression of k-mer sets, WABI 2020 .

[5]: Paolo Ferragina and Giovanni Manzini. "Indexing compressed text." Journal of the ACM (JACM) 52.4 (2005): 552-581.

[6]: Simon Gog, Timo Beller, Alistair Moffat, and Matthias Petri. SDSL-lite: Succinct Data Structure Library, 2016. https://github.com/simongog/sdsl-lite) .

[7]: Yoshihiro Shibuya, Djamal Belazzougui and Gregory Kucherov. locom:Compressing k-mer count tables through minimizers and compressed static functions, 2021. https://github.com/yhhshb/locom .

[8]: Julian Seward. bzip2, 1996. http://www.bzip.org/ .

[9]: Yann Collet. LZ4 - Extremely Fast Compression algorithm, 2011. https://github.com/lz4/lz4 .

[10]: Yann Collet. Zstandard - Fast real-time compression algorithm, 2015. https://github.com/facebook/zstd .

[11]: Armando J. Pinho and Diogo Pratas. MFCompress: a compression tool for FASTA and multi-FASTA data. Bioinformatics, 30(1):117–118, 2013. https://bioinformatics.ua.pt/software/mfcompress/ .

[12]: Shubham Chandak, Kedar Tatwawadi, Idoia Ochoa, Mikel Hernaez, and Tsachy Weissman. SPRING: a next-generation compressor for FASTQ data. Bioinformatics, 35(15):2674–2676, Dec 2018. https://github.com/shubhamchandak94/Spring .

[XXXXXXXXXXXXXXXXX]:

