#Setting up Hestia for the analysis
cd /data/hestia/ksaha
pwd
#----------------------------------------
rm ./CNV_crams/crams/*.cram ./CNV_crams/crams/*.cram.crai    #delete crams that have already been processed
ls -lah ./CNV_crams/crams  #check that all crams have been deleted and there's space for new crams
#----------------------------------------
mv ./CNV_crams/MANTA_ub_pilot/annotation /data/hestia/ksaha/  #move gene annotation to main directory 
rm -r ./CNV_crams/MANTA_ub_pilot #remove the pilot flder to make space
#Check at every step using ls -lah
#-----------------------------------------
mv CNV_crams CNV_references #renames CNV_crams as CNV_references
mv ./CNV_references/crams ./      #move crams to the main directory

#=====================================================================================================================

#Setting up Hestia folders with cram and crai files
#Transfer files using WINSCP to Hestia
ls -lh ./crams #check that all of your crams have been transferred

#LOOP STARTING-------------------------------------------------------------------------------------------------------
for dir in ./crams/CATPD_*; do    
#because all crams and corresponding crai are stored in respective sub-directories, we create a loop to check if anything of the cram or crai is missing
    [ -d "$dir" ] || continue
#-d checks whether the path taken by dir is a directory, if not it skips it
    cram=$(find "$dir" -maxdepth 1 -type f -name "*.cram" -print -quit)
#here cram is a variable that takes the directory's path and searches its files (maxdepth 1 indicates not going to any subdirectories in dir, -type f suggests regular files, and -name *cram finds all anmes ending in .cram
#- print prints the complete path, -quit quits the search after first match is found

    if [ -z "$cram" ]; then
#-z checks hether a string is empty, so if no cram directory was found, then 
        echo "$(basename "$dir"): MISSING CRAM"
#it will echo (output) the name of the directory and missing cram
    elif [ -f "${cram}.crai" ] || [ -f "${cram%.cram}.crai" ]; then
#if a -z suggests cram was found, it will find cram.crai || this stands for or, or it will find file ending with .crai for that sample ID
        echo "$(basename "$dir"): OK"
#if both found, it will say ok
    else
        echo "$(basename "$dir"): MISSING CRAI"
#if not it will say missing crai
    fi
#closes the if statement
done
#closes the loop: the loop goes on until all crams in dir are found
#LOOP ENDING----------------------------------------------------------------------------------------------------------

#=====================================================================================================================

#Using samtools to check if the file structures are okay and can be decoded
#LOOP STARTING---------------------------------------------------------------------------------------------------------
for dir in ./crams/CATPD_*; do
    [ -d "$dir" ] || continue
    cram=$(find "$dir" -maxdepth 1 -type f -name "*.cram" -print -quit)
#same as above where a loop is created for each cram file
    if [ -z "$cram" ]; then
        echo "$(basename "$dir"): MISSING CRAM"
#if no cram, then it will give missing
    elif samtools quickcheck -v "$cram"; then
#otherwise samtools quickcheck will check that cram file's header and whether the expected end-of-file structure is present.
        echo "$(basename "$dir"): QUICKCHECK OK"
    else
        echo "$(basename "$dir"): QUICKCHECK FAILED"
    fi
done
#LOOP ENDING----------------------------------------------------------------------------------------------------------

conda activate manta #activate manta environment
which configManta.py #check which configuration of Manta is running in your active environment
mkdir -p ./results/batch_1 #create results directory 

#======================================================================================================================

#MANTA setup
#Manta is a tool that detects structural variants (SVs) from sequencing data. For our analysis, we want to identify potential deletions and duplications across the whole genome.
#Requirement
#Input sequencing files: The reads from each individual CRAM files, each with a corresponding CRAI index
#Reference genome: The genome against which the reads are aligned	hg38 reference: Homo_sapiens_assembly38.fasta
#Reference index: Allows software to access the reference efficiently: .fai file 
#Manta software: The program that detects structural variants	Manta 1.6.0, installed in the manta Conda environment
#Output directory: Where each sample's results will be stored- 	./results/batch_1/
#Computing resources: CPU resources used to process each sample	-j 8, allowing up to 8 concurrent jobs within a sample's workflow

#define the locations of all files required for MANTA to run
BASE="/data/hestia/ksaha"
CRAM_DIR="$BASE/crams"
OUT="$BASE/results/batch_1"
REF="/home/ksaha/CNV_analysis_manta/Homo_sapiens_assembly38.fasta"

#LOOP STARTING---------------------------------------------------------------------------------------------------------
for dir in "$CRAM_DIR"/CATPD_*; do
#dir is a variable that goes through each sample cram directory as CATPD_* takes all directories starting with CATPD
    [ -d "$dir" ] || continue
    sample=$(basename "$dir")
#basename extracts the sample ID from the sample directory path
    CRAM="$dir/$sample.cram"
#for that directory, it sets CRAM as the cram file in the directory
    RUN="$OUT/$sample"
#$OUT provides the result path, and $sample suggests creating a new directory inside results for that sample, where results will be stored
#----------------------------
#in the fllowing part, first we check whether we miss any input, or we have already had that output before
    if [ ! -f "$CRAM" ]; then
#if there is no CRAM file path for that directory
        echo "$sample: MISSING CRAM — skipping"
        continue
#do not focus on missing crams, continue to next one
    fi
#closes the if clause
#----------------------------
    if [ -e "$RUN" ]; then
#if the output for that sample is already present, 
        echo "$sample: OUTPUT DIRECTORY ALREADY EXISTS — skipping"
        continue
#do not repeat and continue
    fi
#closes the if clause
#----------------------------
#Once these checks are done, if there are no missing inputs and existing output:
    echo "Configuring Manta for $sample"
#interface shows configuring Manta for that sample
#----------------------------
    if ! configManta.py \
        --bam "$CRAM" \
        --referenceFasta "$REF" \
        --runDir "$RUN"; then
        echo "$sample: CONFIGURATION FAILED — skipping"
        continue
    fi
#the above if function is to configure MANTA, configManta.py is a Python script that comes with the MANTA tool, we configure it
#we supply it with the cram file using --bam 
#with the reference fasta file using --referenceFasta
#tell MANTA where it should create the workflow, intermediate files etc using --runDir 
#In this step we have configured Manta by providing everything that it requires
#the ! before specifies that if the configuration fails, it willbe echoes and continued 
#--------------------------
#Once the configuration is done, the interface shows:
    echo "Running Manta for $sample"
#--------------------------
#The following commands are for running the actual MANTA analysis
    if "$RUN/runWorkflow.py" -j 8; then
        echo "$sample: MANTA COMPLETED"
    else
        echo "$sample: MANTA FAILED — check its workflow logs"
    fi
#Here, runWorkflow.py is an executable python script- it reads the workflow configuration set by us, and then runs the workflow
#configManta.py creates the sample-specific workflow directory, and runWorkflow.py is available inside it.
#-j 8 specifies that MANTA can execute a maximum of 8 jobs concurrently. We can change thenumber of concuurent workflows based on computer memory and free usage
#Run following commands for that (no part of the loop)
#nproc
#free -h
#the if clause here checks if the workflow is running for the respective sample
#if yes, it shows manta completed, if not it shows failed
#--------------------------
done
#remember that the entire loop needs to be run in a go, this is where the loop ends
#LOOP ENDING------------------------------------------------------------------------------------------------------------

#=======================================================================================================================

#check if all results there
#LOOP STARTING-----------------------------------------------------------------------------------------------------------
for dir in ./results/batch_1/CATPD_*; do
    [ -d "$dir" ] || continue
    sample=$(basename "$dir")

    if [ -s "$dir/results/variants/diploidSV.vcf.gz" ] &&
       [ -s "$dir/results/variants/diploidSV.vcf.gz.tbi" ] &&
       [ -s "$dir/results/variants/candidateSV.vcf.gz" ] &&
       [ -s "$dir/results/variants/candidateSmallIndels.vcf.gz" ]; then
        echo "$sample: OUTPUT FILES PRESENT"
    else
        echo "$sample: MISSING OR INCOMPLETE OUTPUTS"
    fi
done
#LOOP ENDING-----------------------------------------------------------------------------------------------------------


