# Project summary

Bladder tumors were implanted orthotopically into 3 WT and 4 GSTT2-KO mice (7 samples) at 3-4 months of age and then treated with 4 instillations of M. bovis BCG, following which the bladders were harvested and isolated as single cells for scRNA-seq.

The 2 scRNAseq sublibraries were prepared via the Parse Biosciences Evercode Whole Transcriptome Mini kit, which consists of a total of 12 wells per kit. The sample loading specification table is stored in `sampleLoadingTable.txt`. The output logs from running the various steps of the `Parse biosciences v1.0.4p` pipeline are kept
in the `splitpipe_logs` folder. Old reports submitted to the investigators are 
placed in `reports/`. 

Following the generation of the cells x gene count matrix, the data was analysed in R using the `Seurat` package. Standard downstream Seurat analyses was conducted, including dimensional reduction and clustering, marker gene identification, manual cell type annotation and pathway enrichment analysis. Scripts for differential abundance (DA) analysis and differential state analysis are included alongside the aforementioned analysis in `scripts/`, which were run in a previous round of analysis. 

# Location of raw data

s3://parse-biosciences-mugdha
On Biodebian, the files are stored at: /media/gedac/kane/projects/parse-bladder-scrnaseq

# Download

The files in this repository are also stored in an S3 bucket: s3://parse-scrnaseq-bladder/. The raw and processed data files are kept in the 
`data/` folder. R objects are kept in the `rds/` folder. 

# Contact
PhD student: [Mugdha Vijay Patwardhan](mugdha.p@u.nus.edu)
Principal Investigator: [Ratha Mahendran](surrm@nus.edu.sg)
Bioinformatician: [Toh Qin Kane](kane9530@hotmail.com)

# Chargeable hours

Approximately 15h (extensive amount of work, approx 80+h since October 2022,was performed as a gesture of goodwill, prior to an official consulting model being put in place).

