# A.C.Rsuite
ATAC-, ChIP- and RNA-seq data analysis tool suite.

This is R (a programing language) based command line tool (CLT) suite which provides the preliminary graphic results as well as spreadsheets. It also tries to produce the publication-ready figures.
It is mostly for wetlab scientists to generate the preliminary results, and to get a scene if the experiments worked before handling to bioinformatician.
There is no install step, download or clone, as well as put the path in the system "PATH", and then enjoy!

And the wetlab user can also learn data processing using R through all the scripts used here.  

[Pipeline Figure coming soon...]

The CLT was initially built in Chris Class lab @UCSD. If you have questions and requests, please submit an issues. Enjoy!
# Pre-requirements
  - R,
  - [***homer***](http://homer.ucsd.edu)

# Instruction
There are two main pipelines, one is for standard RNA-Seq analysis; the other one is for ATAC-Seq and transcription factor (TF) ChIP-Seq. Currenly it dose not support broad peak assays such as histone ChIP. But you can use peaks from ATAC-seq to quantify them, such as H3K27ac ChIP-Seq.

All commands below provide a help instruction by no parameter, such as `rnaPipe.R` in command line.

## Prepare the sample definition file
In order for the pipelines to know the locations of your samples of all groups/conditions, the sample definition file is required. It is a 4 column table for RNA-Seq or ATAC-Seq and 5 column table for ChIP-Seq. The table is 'tab' (*\t*) separated without header. Each row is a group.

- First column: the name of the group
- Second column: the color of the group in the plots, [R-color](http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf) by name can be used, as well as [RGB](https://colorbrewer2.org/) code as *#de2d26*.
- Third column: the full path to the tag directories of samples/replicates related to this group, separated by semicolon (*;*)
- Fourth column: the name of each sample, separated by semicolon (*;*) as well, the number of sample names should be the same as the number of tag directories

For TF ChIP-seq, there is one addition column:
- Fifth column: the full path to the input tag directories corresponding to the samples in column 3, separated by semicolon (*;*), the number should be the same as the number of tag directories. You can repeat the input tag directory if needed

## rnaPipe.R
```
rnaPipe.R /path/to/sample/definition/file [-o /path/to/output/folder/ -g /genome/ -l /gene/length/cutoff -t /UCSC/genome/track/name -m /minimal/TPM/cutoff -c /additional/homer/commands/for/annotation/calling]
```
*[]* means the optional parameters, you can choose which one/ones to use, but don't include the *[]* in your command.
- /path/to/sample/definition/file, is required, and the format is introduced above
- /path/to/output/folder/, the full path to the output directory
- /genome/, hg38 or mm10, which has to be installed in ***homer***
- /gene/length/cutoff, the numeric length cutoff for the small genes which will not be considered.
- /UCSC/genome/track/name, the UCSC genome browser track will be generated, and can be added as "/homer_data/www/html/hubs/***TRACK NAME***/trackDb.txt"
- /minimal/TPM/cutoff, the numeric minimal TPM cutoff,

### RNA-Seq Pipeline Steps
The RNA process pipeline includes the following steps. The user can choose any step to run or rerun. Some steps contain more options to use. You can type the command of each step, and a usage help with available options will show in the command line.

1. `alignStats.R`, produces *alignStats.txt* file in the output directory. It contains all the mapping QC. This helps you to QC the sequencing.
2. `rnaQuan.R`, produces raw count matrix and matching TPM matrix by ***homer*** in the folder named *rnaQuan* in the output directory.
    - Obtain the raw count matrix by '-condenseGenes -count exons -noadj' options with ***homer***;
    - Obtain the TPM without '-condenseGenes' options, then subset it with genes reported by above raw count matrix;
    - Remove the short genes, as well as plot the gene length agaist mean and sd of expression, and *rawT.txt* as well as *rawC.txt* will be saved in *rnaQuan* folder;
    - Plot the principle component analysis (PCA) of all samples with specified grouping colors defined in sample definition file in *overall.PCA.pdf*;
    - Calculate and plot the pair-wised TPM expression Pearson's correlation of all samples, as well as pair-wised TPM expression Pearson's correlation among replicates of each group. This is very important and useful to check the quality of the replicates;
    - Merge all replicates of a group to generate a tag directory per group in *mergeTag* folder;
    - Make a UCSC genome browser hub file on homer server, user can easily add hub on UCSC by <div style="display: inline">http://homer.ucsd.edu/hubs/HUB NAME/hub.txt</div>
3. `rnaDiff.R`, perform the differential expressed gene (DEG) analysis by DEG. It performs group pair-wised, **NOT** all groups together. The reason is that the pair-wised results are independent with how many groups are included as well as sample changes in other groups. And the results should NOT be much difference.
    - Additional options, such as logFC cutoff, FDR cutoff, are available, and **NO** need to run the whole pipeline again;
    - Please provide the *rawC.txt* and *rawT.txt* from *rnaQuan* folder to the DEG analysis.














The end
