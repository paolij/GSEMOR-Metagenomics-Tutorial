# UF 2025 GSEMOR Workshop
# Metagenomic Classification Tutorial
# Julia Paoli

### **Background**

In this section on metagenomics for pathogen identification, we will investigate the viral and bacterial species present in sequencing data. This is an important initial step in an outbreak scenario to identify and characterize the pathogens present in your samples. 

The data we will analyze comes from a surveillance study of ticks collected in Mongolia. We collected *I. persulcatus* ticks (n = 1300) in May 2020. Viral RNA was extracted from cell cultures, converted to cDNA, and library prepped. We conducted whole-genome shotgun sequencing on the Illumina iSeq100 platform. To read more about this study, see the paper [here](https://www.mdpi.com/2076-0817/13/12/1086).

We will be submitting jobs to the HPG server using SLURM job scripts. For more information on submitting jobs check out this [website](https://help.rc.ufl.edu/doc/HPG_Scheduling). For good computing practices using HPG check [here](https://help.rc.ufl.edu/doc/HPG_Computation). 

⚠️ Warning: Ensure you are not running analyses on login nodes. To avoid this, submit jobs to the scheduler using sbatch command followed by your script.

### **Kraken Metagenomic Pipeline**

First, we will use the Kraken suite to classify, quantify, and visualize our metagenomic dataset following the steps outlined in [this paper](https://www.nature.com/articles/s41596-022-00738-y).

1. Taxonomic classification: Kraken2
2. Abundance Calculation: Bracken
3. Microbiome visualization: KrakenTools

Let's get started!

1. **Connect to HiPerGator (HPG) server to run analyses**

    Navigate to your preferred terminal and enter the following code:

    ```bash
    ssh <GatorLink Username>@hpg.rc.ufl.edu
    ```

    When prompted, enter your Gatorlink password and duomobile authentication

2. **Set up your folders**
    
    Navigate to the workshop directory in the blue server

    ```bash
    cd /blue/general_workshop
    ```

    Navigate to your folder in the workshop directory using the command cd. If you haven't made one yet, please take the opportunity to do so now using the command mkdir <your_name>*

   ```bash
   cd <your_folder>
   ```

    Make a folder for metagenomics and go into newly created folder
    ```bash
    mkdir metagenomics
    ```
    ```bash
    cd metagenomics
    ```

4. **Run Kraken2 to classify metagenomic reads**


   HPG already has Kraken databases available for use so no need to download or build them.

   Let's write a script!

   ⚠️ Please remember to update the email parameter in the script to reflect your own information, I've put in place holders for now.
   
   ⚠️ When running paired reads with Kraken2, make sure to keep the # symbol in the classified-out option for proper handling of paired data
   
    ```bash
    nano kraken2.sh
    ```

    You can copy this script and paste it to your newly opened file.
    To save your script do Ctrl + O then Ctrl + X to exit.
    ```bash
    #!/bin/bash
    #SBATCH --job-name=KRAKEN2
    #SBATCH --mem=10G
    #SBATCH --time=24:00:00
    #SBATCH --nodes=1
    #SBATCH --cpus-per-task=1
    #SBATCH --account=general_workshop
    #SBATCH --qos=general_workshop
    #SBATCH --mail-user=<username>@ufl.edu #put in your email
    #SBATCH --mail-type=ALL
    
    pwd; hostname; date

    # Create output directory
    mkdir -p kraken_output

    # Load Kraken2 module
    module load kraken/2.1.3
    
    # Input arguments for paired-end reads
    FWD=$1
    REV=$2
    
    # Directory for Kraken2 database
    KRAKEN_DB_PATH="/data/reference/kraken2/minikraken2_v2_8GB_201904_UPDATE"
    
    # Extract the base name of the forward reads file without path and extension
    BASE_NAME=$(basename "$FWD" | sed 's/\..*//')
    
    # Run Kraken2
    kraken2 --db $KRAKEN_DB_PATH \
            --paired $FWD $REV \
            --classified-out kraken_output/${BASE_NAME}_classified_seqs#.fq \
            --output kraken_output/${BASE_NAME}_kraken_output.out \
            --report kraken_output/${BASE_NAME}_kraken_report.out \
            --threads 8 \
            --minimum-hit-groups 3
    ```

    To submit the job use sbatch. 

    In the script, we assigned the forward read to the value $1 and the reverse read to the value $2. These are positional parameters that specify the order of files you provide as     command-line arguments. The general format for submitting a job is:
    
    sbatch script.sh $argument1 $argument2

    We will use the trimmed sequence read files in the shared workshop folder as input for Kraken2. We can indicate the path of files for input in the script:

    ```bash
    sbatch kraken2.sh /blue/general_workshop/share/metagenomics/6_S6_L001.fwd_p.fq.gz /blue/general_workshop/share/metagenomics/6_S6_L001.rv_p.fq.gz
    ```
    To check on the status of your job
    ```bash
    squeuemine
    ```

    If you look at the slurm output file you will see a percentage breakdown of the classified and unclassified reads.

6. **Run Bracken for abundance estimation**

    We will now use the output from Kraken2 as input for Bracken abundance estimation.

    ```bash
    nano bracken.sh
    ```

    ```bash
    #!/bin/bash
    #SBATCH --job-name=Bracken
    #SBATCH --mem=10G
    #SBATCH --time=24:00:00
    #SBATCH --nodes=1
    #SBATCH --cpus-per-task=1
    #SBATCH --account=general_workshop
    #SBATCH --qos=general_workshop
    #SBATCH --mail-user=<your username>@ufl.edu #put in your email
    #SBATCH --mail-type=ALL

    pwd; hostname; date

    # Create output directory
    mkdir -p bracken_output

    # Load bracken module
    module load bracken/2.9

    # Directory for Kraken2 database
    KRAKEN_DB_PATH="/data/reference/kraken2/minikraken2_v2_8GB_201904_UPDATE"

    # Sample Input Argument
    FWD=$1
    
     # Extract the base name of the forward reads file (without path and extension)
    BASE_NAME=$(basename "$FWD" | sed 's/\..*//')

    # Run Bracken
    bracken -d $KRAKEN_DB_PATH -i kraken_output/${BASE_NAME}_kraken_report.out -r 100 -l S -t 10 -o bracken_output/${BASE_NAME}.bracken -w bracken_output/${BASE_NAME}.breport
    ```

    To run the script pass the forward read as the first argument. We are doing this to extract the basename of our sample. 
    ```bash
    sbatch bracken.sh /blue/general_workshop/share/metagenomics/6_S6_L001.fwd_p.fq.gz
    ```

7. **Visualize Microbiome Classification and Abundance**

    We will now make Krona plots using the output from Bracken to visualize the microbiome composition of our data in pie chart format.
    
    ```bash
    nano krona_plot.sh
    ```

    ```bash
    #!/bin/bash
    #SBATCH --job-name=Krona_Plot
    #SBATCH --mem=1G
    #SBATCH --time=24:00:00
    #SBATCH --nodes=1
    #SBATCH --cpus-per-task=1
    #SBATCH --account=general_workshop
    #SBATCH --qos=general_workshop
    #SBATCH --mail-user=<your username>@ufl.edu #change your email
    #SBATCH --mail-type=ALL
    
    # Load required modules
    module load krona/20161019
    module load krakentools
    
    # Sample Input Argument
    FWD=$1

    # Extract the base name of the forward reads file (without path and extension)
    BASE_NAME=$(basename "$FWD" | sed 's/\..*//')
    
    # Create output directories
    mkdir -p b_krona_txt
    mkdir -p krona_html
    
    # Generate Krona input from bracken report
    python /blue/general_workshop/share/metagenomics/kreport2krona.py -r bracken_output/${BASE_NAME}.breport -o b_krona_txt/${BASE_NAME}.b.krona.txt --no-intermediate-ranks
    
    # Create Krona HTML plot
    ktImportText b_krona_txt/${BASE_NAME}.b.krona.txt -o krona_html/${BASE_NAME}.krona.html
    ```
    
     To run the script pass the forward read as the first argument. We are doing this to extract the basename of our sample. 
    ```bash
    sbatch krona_plot.sh /blue/general_workshop/share/metagenomics/6_S6_L001.fwd_p.fq.gz
    ```
    
    To view the plot, download the html file using your preffered file transfer system and open in your internet browser.
   

**BLAST Metagenomic Pipeline**

BLAST (Basic Local Alignment Search Tool) is a powerful tool to find regions of similarity between nucleotide or protein sequences. The user inputs a "query" sequence and the program searches a large database of nucleotide or protein sequences to find the best match to the query sequence. BLAST is useful for identifying your study contigs and finding similar sequences. 

While BLAST has an easy to use web [interface](https://blast.ncbi.nlm.nih.gov/Blast.cgi), running BLAST through the command line is easier when dealing with a large query set to automate the process. The command lines also allows for customization of the search database and is often faster than running on the web interface. 


1.**Write BLAST Script**

Let's write a script to compare our contigs to the NCBI viral database and find the best match. We will use BLAST output [format 6](https://www.metagenomics.wiki/tools/blast/blastn-output-format-6) to provide information to assess how good each hit is to our query.

```bash
nano blast_virus.sh
```

```bash
#!/bin/bash
#SBATCH --job-name=BLAST_VIRUS
#SBATCH --mem=10G
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --account=general_workshop
#SBATCH --qos=general_workshop
#SBATCH --mail-user=<your username>@ufl.edu
#SBATCH --mail-type=END,FAIL

# Load BLAST module
module load ncbi_blast

# Create output directory if it doesn't exist
mkdir -p blast_virus

# Input FASTA file
FILE="$1"
BASENAME=$(basename "$FILE" .fasta)

# Run BLASTn
blastn -db ref_viruses_rep_genomes \
       -query "$FILE" \
       -out "blast_virus/${BASENAME}.txt" \
       -num_threads 8 \
       -outfmt "6 qseqid sacc stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"
```

To run the script, we will use as input a final contig file generated by de novo assembly using MEGAHIT. We will pass the fasta file as the first argument.

```bash
sbatch blast_virus.sh /blue/general_workshop/share/metagenomics/6_S6_L001_final.contigs.fa
```

Now you can download the text file generated to look at your BLAST hits. 

To assess quality of the BLAST hits you can assess:

1. ***Evalue***: the closer to 0 the better the quality of the match. BLAST E-value is the number of expected hits of similar quality that could be found just by chance.

2. ***Percent Identity***: the closer the value to 100% the more similar your query sequence is to the BLAST hit. Calculated by dividing the number of identical nucleotide matches in an alignment by the length of the alignment.

3. ***Bit-score***: the higher the bit-score the better the more similar the two sequences are. Bit-scores are normalized and so independent of database search size.

For a nice glossary of BLAST terminology see this [glossary](https://www.ncbi.nlm.nih.gov/books/NBK62051/)
