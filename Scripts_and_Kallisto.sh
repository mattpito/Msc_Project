Script to iterate fastqc. Make sure the module fastqc is loaded in order for it to work

#!/bin/bash
#This code enters the directory nm001 in the birkbeck server where the downloads
#are located and is using them as input. We target the files with .gz extension
echo $'Now running fastqc'
#the code executes in a sub-shell. Will enter the directory below and iterate.
#CHANGE THE PATH TO WHERE YOUR FILES ARE LOCATED
cd /d/projects/u/nm001/brain_tumour/acp_samples/
for i in *.gz;do
    echo 'Running FastQC with '$i
    fastqc $i -o /d/projects/u/pm005/msc_project_2021/michele_fastqc/
done

#######################################################################################################
Script to isolate read counts from our files. 
 
 
#!/bin/bash

cd /d/projects/u/nm001/brain_tumour/acp_samples/
for i in *.gz;do
    printf "%-30s | %-30s \n" "$i" "$(($(zcat $i| wc -l) / 4))" >> /d/projects/u/pm005/msc_project_2021/read_counts.txt

  
done
#######################################################################################################

Code to compare reads (previous scrip output) with the ENA database in R :
# open read_counts as dataset x
x = read.table(file ='FILE_PATHWAY', header = FALSE, sep = "|", dec = ".")
# downloaded TSV file of the ENA read counts and stored them as y dataset
y = read.table(file ='FILE_PATHWAY', header = TRUE, sep = "\t", dec = ".")
#The paired entries from ENA are listed as one. I duplicate them to match each reading individually. Finally I add the ENA counts to x dataset.
y2 <- y[rep(seq_len(nrow(y)), each = 2), ]
y3<-subset(y2,select=c(read_count))
x <- cbind(x, userID = y3$read_count)
library(knitr)
kable(x,align='lll',
col.names = c("Entries","Read_Counts",'ENA_Read_Counts'))

#########################################################################################################
replace FILE with your version. Genarates index on the reference genome as well as output path
kallisto index -i gencode.v37.transcripts.fa.index /d/projects/u/pm005/msc_project_2021/reference_gen/FILE.fa.gz


This script runs kallisto on all files to generate abundancies. Replace the FILE with your generated index and output path
#!/bin/bash
cd /d/projects/u/nm001/brain_tumour/acp_samples/
for file in *_1.fastq.gz;do
    name=$(echo "$file"| sed 's/_.*//' )
    echo "Now working on paired entry" $name
    kallisto quant --rf-stranded --threads=8 -b 100 --index=/d/projects/u/pm005/msc_project_2021/reference_gen/FILE.fa.index --output-dir=/d/projects/u/pm005/msc_project_2021/kallisto/$name ${name}_1.fastq.gz ${name}_2.fastq.gz
done









