Welcome to FC-Viurs, an algorithmic software for full-length genome assembly for viruses.
The package contains binary files that can be run directly.

Requirements
Linux operating system
c++11 or higher

Usage
1. Download FC-Viurs
www.github.con/FC-Virus download zip file or download directly
Clone the repository:
   ```bash
    "git clone https://github.com/TheyuGao/FC-Virus.git"
2. Unzip the FC-Virus package (make sure you have unzip installed)
Execute the command in the FC-Virus download directory
    "unzip -l FC-Virus-master.zip"
3.install FC-Virus
Generally speaking, FC-Virus can be installed using the Make command.
First, go to the FC-Virus directory:
    "cd /path/in/FC-Virus"
Then install
    "make"
4.Add environment variables
 if you want to run FC-Virus in any directory, add environment variables.
    "nano ~/.bashrc"
    "export PATH="$PATH:/path/to/your/project/bin""
    "source ~/.bashrc"
   
5. Run FC-Virus (requires std=c++11)
If you want to skip or run FC-Virus directly after completing the above step4
"cd FC-Virus/bin"
"FC-Virus" or "./FC-Virus"
For novice users, we recommend you to run our binary file "FC_Virus" directly.
in the directory /your PATH/FC-Virus/bin
". /FC-Viurs -o PATH(your/output/path) -k k_number(default 25) --left PATH(/your/fq_1/file) --right PATH(/your/fq_2/file) -fr fr_number(1, 2or3) "

For those of you who have other needs, we recommend that you compile and run it directly,The source code is stored in /FC_Virus/code.

g++ -o FC-Virus -std=c++11 main.cpp GeneralSet.cpp kmer.cpp ex_r.cpp ex_l.cpp

Example
". /FC-Viurs -o test -k 25 --left . /test/fq1.fastq --right . /test/fq2.fastq --fr2"

The contig results are output in the specified -o output path, named "FC_contig.fa".
