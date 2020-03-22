# Pipelines and scripts


This scripts were made to analyse the proportion of new and old genes in different phases of spermatogenesis. The scripts were developed during my mastersd research project from 2015 to 2017. 
The [final.pl](final.pl) perl script needs three inputs: 

- the supplemental table 1 from Vibranovski _et al_ (2009) in PLoS Genetics. This table has the expression of genes in different spermatogenesis' stages. The table has 3 replicates for each stage of spermatogenesis;

- two suplemental tables from Zhang _et al_ (2010) in Genome Research. The first table has the age of _Drosophila_ genes, as the branch in which they appear in the phylogeny (supplemental table S1). The second table from this paper has values of dN and dS for such genes (supplemental table S4a).

The [final.R](final.R) and [positiveselection.R](positiveselection.R) R scripts need the same two inputs, and they do not have to be inputed when calling the scripts:

- the output from the [final.pl](final.pl) script, called final.output by default;

- [flyDIVaS](http://www.flydivas.info/) data for the melanogaster subgroup, called melsubgroup_analysis_results_flydivas, version 1.2.

The last R script, [simulation.R](simulation.R) has two functions simulating the time in generations needed for a new allele to be fixed in a haploid and diploid populations. The difference betweeen the functions is that one allows for genetic drift, while the other one is deterministic. The code for both was modified from http://rosetta.ahmedmoustafa.io/selection/ . The script doesn't need any input file, and it automatically runs with pre-set values that can be changed in the last lines of the script.
Finally, the bash script [images.sh](images.sh) just converts the pdf figures created by the R scripts into png and jpg images.


## Getting Started

To get the  scripts running you need to just download the scripts and run them according to the usage instructions. You can get the instructions for each script by reading bellow, or by running them without any inputs or parameters.

### Prerequisites

To be able to run the scripts you need bash, Perl version 5.18 and R version 3.2. Please follow the developers' directions for download and installation of the programs.

R: https://cloud.r-project.org/

Perl: https://www.perl.org/get.html


## Usage

For using the *final.pl* with an input with the age data, one with _dN_ and _dS_ data, and one input with spermatogenesis expression:

`
perl final.pl -exp EXPRESSION_TABLE -dnds DNDS_TABLE -age AGE_TABLE
`

The files used as input and their origins can be found in the above session "Pipelines and scripts".
After running the script, two files will be created: a log file, and an output file. The log is called *final.log*, and the output is *final.txt*. The log file will have a compilation of all the runs for the script, and any errors or problems with the run.

Once the output file is produced it is possible to run both [final.R](final.R) and [positiveselection.R](positiveselection.R):

`
Rscript final.R
`

`
Rscript positiveselection.R
`

The [simulation.R](simulation.R) can be ran at any point since it takes no input files. It is still possibleall of that is done, it is possible to run the R script to get statiscal analysis and plots.

`
Rscript simulation.R
`

Or run the following on R:

`
source("simulation.R")
`

`
diploid_haploid_deriva(selection_coeficient, number_of_generations, number_of_simulations, frequency_recessive, population_size) # with drift
`

`
diploid_haploid_determinista(selection_coeficient, number_of_generations, frequency_recessive) # deterministic
`

Once all the scripts have been ran, it is possible to run [images.sh](images.sh) to get .png and .jpg versions of the .pdf images.


`
./images.sh
`


## License

This project is licensed under copyleft agreement. You may use it and modify it but should cite this project and author when publishing.

## Acknowledgments

* Prof. Dr. Maria Vibranovski (advisor)
