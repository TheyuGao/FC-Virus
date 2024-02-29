Welcome to FC-Viurs, an algorithmic software for full-length genome assembly for viruses.
The package contains binary files that can be run directly.

Requirements
Linux operating system
c++11 or higher

Usage
1. Download FC-Viurs
www.github.con/FC-Virus download zip file
2. Unzip the FC-Virus package (make sure you have unzip installed)
Execute the command in the FC-Virus download directory
"unzip -l FC-Virus-master.zip".
3. Run FC-Virus (requires std=c++11)

For novice users, we recommend you to run our binary file "FC_Virus" directly.
in the directory /your PATH/FC-Virus/bin
". /FC_Viurs -o PATH(your/output/path) -k k_number(default 25) --left PATH(/your/fq_1/file) --right PATH(/your/fq_2/file) -fr fr_number(1, 2or3) "

For those of you who have other needs, we recommend that you compile and run it directly,The source code is stored in /FC_Virus/code.

g++ -o filename -std=c++11 main.cpp GeneralSet.cpp kmer.cpp ex_r.cpp ex_l.cpp

Example
". /FC_Viurs -o test -k 25 --left . /test/fq1.fastq --right . /test/fq2.fastq --fr2"
". /filename -o test -k 25 --left . /test/fq1.fastq --right . /test/fq2.fastq --fr2"

The contig results are output in the specified -o output path, named "FC_contig.fa".
