# FC-Virus

Welcome to **FC-Virus**!  
FC-Virus is a powerful tool designed for full-length virus genome assembly. The package includes precompiled binaries for direct execution.

### Requirements:
- Operating System: Linux
- C++ Standard: C++11 or higher

### Installation Instructions:
We have submitted the bioconda go-live application, so stay tuned.
After the application is approved you can install it with the following command:
1. **Install via Bioconda**:  
   Simply run the following command:
    ```bash
   bioconda install FC-Virus
     ```

3. **Use Precompiled Binaries**:  
   Download FC-Virus from GitHub or clone the repository:
    ```bash
   git clone https://github.com/TheyuGao/FC-Virus.git
    ```
   Ensure you have unzip installed. Unzip the downloaded file:
   ```bash
   unzip FC-Virus-master.zip
   ```
   Navigate to the FC-Virus directory:
   ```bash 
   cd /path/to/FC-Virus
   ```
   Install by running:
   ```bash 
   `make`
   ```
   Set Up Environment Variables: To run FC-Virus from any directory, add it to your PATH:
   ```bash 
   nano ~/.bashrc
   ```
   Add the following line:
   ```bash 
   export PATH="$PATH:/path/to/FC-Virus/bin"
   ```
   Apply the changes with:
   ```bash 
   source ~/.bashrc
   ```

5. **Compile from Source**:  
   To compile from source, navigate to the code directory:
   ```bash 
   cd /path/to/FC-Virus/code
   ```
   Compile using g++:
   ```bash 
   g++ -o FC-Virus -std=c++11 main.cpp GeneralSet.cpp kmer.cpp ex_r.cpp ex_l.cpp
   ```
   Then run FC-Virus:
   ```bash 
   ./FC-Virus -o /path/to/output -k 25 --left /path/to/left.fq --right /path/to/right.fq
   ```

### Example Command:
```bash 
./FC-Virus -o test -k 25 --left /path/to/test/fq1.fastq --right /path/to/test/fq2.fastq
```

The contig results will be saved to the specified output path as `FC_contig.fa`.
