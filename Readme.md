CancerCoCoPUTs V1.1
Author: Douglas Meyer, US Food and Drug Administration

Prerequisites:

Python 3 (3.6 or higher recommended)
Pandas   (1.1.5 or higher recommended)

Included Data Files:

categories.tsv
	provides relevant metadata for each read count file

ensembl_transcript_lengths_w_ENSG.tsv
	primary transcript ENST label and transcript length for
	each of ~19000 human genes (ENSG)

Supplemental Table 1.xlsx
	details which tissue samples are grouped together in 
        each tumor type and normal tissue type described in 
        <paper>
	
GRCh38_<cocop>_w_ENSG.tsv
	4 different files (codons, bicodons, dinucleotides, 
        junction_dinucleotides)
	primary transcript ENST label and the number of <cocop>
        for each of ~19000 human genes (ENSG)
		for example, codons file contains information 
                about the number of each codon present in the 
                primary transcript of each of ~19000 human genes
	
	***NOTE***
	Due to its size, <GRCh38_bicodons_w_ENSG.tsv> is split into
	10 separate files (GRCh38_bicodons_w_ENSG_x.tsv for x=1-10).
	<GRCh38_bicodons_w_ENSG.tsv> can be generated by concatenating
	the 10 files together: 
	cat file1 file2 ... file10 > GRCh38_bicodons_w_ENSG.tsv
		

Other Necessary Files:

READ_COUNTS
	can be downloaded from 
	https://portal.gdc.cancer.gov/repository
		"Workflow Type" = "HTSeq - Counts"



Pipe.py includes all methods necessary to compute transcripts per 
million (TPM) based on one or more RNAseq read count file. It also 
includes all methods necessary to compute transcriptomic weighted 
codon usage, codon pair usage, dinucleotide usage and junction 
dinucleotide usage.

<paper> describes how this code was used to analyze different cancer 
and normal tissue types from RNAseq data available through TCGA and 
how it was used to generate data for the novel database CancerCoCoPUTs

Driver_code.py shows how pipe.py can be used to generate 5 different 
sets of files based on TCGA RNAseq read count data. ?folder? parameter 
(path to folder containing rnaseq read count files) and ?out_dir? 
parameter (path to desired output folder) must be set prior to run. 
Successful run will result in 25 files saved to out_dir, 
including codon usage, codon pair usage, dinucleotide usage, 
junction dinucleotide usage and TPM files for each set type.

5 types of file sets
       All_cancer: each individual tumor sample
       All_normal: each individual normal tissue sample
           Paired: each pair of tumor and normal tissue samples 
                   collected from the same patient
    Median_cancer: values computed based on median TPM for each 
                   tumor type 
    Median_normal: values computed based on median TPM for each 
                   normal tissue type
