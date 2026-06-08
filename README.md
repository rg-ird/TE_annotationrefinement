# TE Annotation Refinement (TEAR)  

Pipeline for plant transposable element annotation refinement.  
This Pipeline use HiTE results (repeatmasker .out and library files), Tesorter results (with the HiTE library files) and the Inpactor 2 results.  
It processes the HiTE results repeatmasker .out file to remove overlapping annotation based on repeatmasker score and annotate the LTR retrotransposons at the family level.  
**Prerequisite:**     
-[HiTE] (https://github.com/CSU-KangHu/HiTE) annotation results: .out file   
-[TEsorter] (https://github.com/zhangrengang/TEsorter) Hite_library annotation results: confident_TE.cons.fa.rexdb-plant.cls.tsv   
-[Inpactor2] (https://github.com/simonorozcoarias/Inpactor2) Inpactor 2 library results: Inpactor2_library.fasta  

# Workflow
 
![Workflow](docs/workflow.png

# Installation 
git clone https://github.com/rg-ird/TE_annotationrefinement.git  
cd TE_annotationrefinement  
You need to install/load NCBI blast before to run the pipeline.  
You also need to sort the chromosome lenth into a file for graphical outputs:  
samtools faidx genome.fa   

# Usage
python tear_pipeline.py --param tear_pipeline_config_template.txt  

# Results  

In TEAR_RESULTS atre the results of each step.  
-the final filtered .out  
-the .tbl file  
-All graphics in PDF format, including the csv table of density.






  

        
