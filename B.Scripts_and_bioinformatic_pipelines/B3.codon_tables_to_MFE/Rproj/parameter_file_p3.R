###############################################################################
###################################BAR 18.01.24################################
########## parameter file for running the VirBio scripts for Relative Enrichment
###############################################################################
# before running, make sure your directory contains:

# CUSTOM FILES
# comparisons_Xdms.csv #custom
# Xdms_library_info.xslx #custom

################################################################################
## select region to run and add path to DMS_analysis_VirBio from your working 
## directory. Then run!
## (make sure to adjust MFE parameters (min.cov, pseudocount) to your analysis 
## below)

path.to.scripts="../"

################################################################################

## for 1. scripts virvar codon table###########

#where your virvarseq codon files live
codon.dir=paste0("../codontables/",region)

# length of the fragment
aa.length=ifelse(region == "p1", 851, 
          ifelse(region == "p2", 628, 
          ifelse(region == "p3", 770, NA)))

# where the fasta file 
fasta.reference.path=ifelse(region == "p1", paste0(path.to.scripts,"/generic_files/p1_dms.fa"), 
                    ifelse(region == "p2", paste0(path.to.scripts,"/generic_files/p2_dms.fa"), 
                    ifelse(region == "p3", paste0(path.to.scripts,"/generic_files/p3_dms.fa"), 
                                                  NA)))

# where the name conversion file is
# has 2 columns:name == file name & ref.name == new name
# name.conversion.file = paste0("./",region,"dms_library_info.xlsx")
name.conversion.file = NA

# where the codon identity file lives
codon.identity.file=paste0(path.to.scripts,"/generic_files/codon_identity_matrix.csv")

# if you introduce specific mutations and not all
# where the file with the expected mutations is

mutations.introduced=ifelse(region == "p1", NA, 
                    ifelse(region == "p2", paste0(path.to.scripts,"/generic_files/mutations_introduced_in_P2.csv"), 
                    ifelse(region == "p3",  paste0(path.to.scripts,"/generic_files/mutations_introduced_in_P3.csv"), 
                    NA)))

#type of analysis and working directory 
working.dir=ifelse(region == "p1", "filter_single_mutants", 
            ifelse(region == "p2", "filter_expected", 
            ifelse(region == "p3", "filter_expected", NA)))

# minimum coverage to call a deletion
min.coverage.deletions=1

# the strand bias cutoff 
sor.cutoff=4

# Filter positions in the table because they are not mutagenized?
filter.positions=T

# AA where mutagenesis start and ends
AA.start.mut=ifelse(region == "p1", 1, 
                     ifelse(region == "p2", 18, 
                     ifelse(region == "p3", 1, NA)))

AA.end.mut=ifelse(region == "p1", 851, 
           ifelse(region == "p2", 605, 
           ifelse(region == "p3", 770, NA)))

# do we correct the position of the fragment to match genomic AA
correct.position=T

# if so, what is the first AA of the mutagenized fragment
startAApos=ifelse(region == "p1", 1, 
                  ifelse(region == "p2", 847, 
                         ifelse(region == "p3", 1416, NA)))

# the min number of reads to keep a mutation in a table
min_reads=0

################################################################################
################################################################################
################################################################################
## parameters for MFE

## file giving the preselection/control and the post/selection
## needs columns of pre, post, name
#add groupings as desired (per prelicate, per lib, per region...)
comparison.file=paste0("./comparisons_",region,"dms.csv")

# what is the pseudo count
pseudo.count=1

### what to filter in pre, post or both 

## you can remove sites that don't get a minimum value in the pre
## value > 0 is good for MFE
## for diff sel could increase accuracy when dealing with negative values but reduce resolution in positive
min.cov.pre=5

## you can remove mutations that don't reach a minimum in pre OR post
#  this can be used in diff selection to make sure you only call mutations that go from a real count to 
# something lower or higher
min.cov.pre.or.post=NA

#minimum coverage for stops
min.cov.pre.or.post.stops=NA

## you can remove mutations that don't reach a minimum at both pre and post
# this can increase accuracy in diff selection analyses
min.cov.pre.post=NA

##  Do we remove stops from the analysis?
rm.stops=T

################################################################################
################################################################################
##parameters for diffsel

### if you want to plot the mean/median cutoff
# default=mean +2*sd
mean.or.median.cutoff="mean"
sd.times.cutoff=2

#cvb3 reference file for protein positions
ref.file= paste0(path.to.scripts,"/generic_files/cvb3_reference_sites.csv")