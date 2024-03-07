Data and scripts for "Mapping the mutational landscape of a full viral proteome reveals distinct profiles of mutation tolerability"


A. Data 

A1.NGS_library_prep_protocol: contains the protocol for the preparation of duplex sequencing adapters (Protocol_duplex_adapters.docx), preparation of NGS dual index libraries (Protocol_duplex_NGS_library_prep.docx) and sequences of primers used (Primers_duplex.xslx) and template for the qPCR standard curve (Luciferase_stnd_curve.dna). 

A2.synthetic_oligonucleotides_P2_P3: sequences of the synthetic oligonucleotides used for the P2 and P3 DMS. 

A3.codon_tables: codon counts of regions P1, P2 and P3 in the plasmid libraries, and passage 2 virus obtained from HeLa-H1 and RPE cells.  

A4.comp_assay_data: data from the competition assays. The ratio of fluorescent counts at 20hpi vs 8hpi of mutant or WT mCherry-CVB3 versus WT GFP-CVB3 virus are indicated. 
  
A5.site_and_mut_tables: contains data at the mutation and site level for the full proteome. Includes site and mutation MFE, DDG, Shannon entropy in enterovirus B alignments, RSA, secondary structure, and annotations. 
  
A6.seq_alignments_and_phydms: aligned and trimmed enterovirus B sequences used for the calculation of Shannon’s entropy, and phydms results per region. 

A7.pymol_sessions_map_MFE: Pymol sessions of high-resolution or Alphafold2-predicted structures colored according to average site MFE values by their assignation to b-factors. 

A8.RPE_vs_Hela_diffsel: Site differential selection values per site in Hela-H1 and RPE. 
 
A9.sitemap_data: Schrodinger SiteMap software output tables for predicted druggable pockets in CVB3 protein structures. 

 

B. Scripts_and_bioinformatic_pipelines 

B1.DMS_oligo_design_pipeline: Pipeline for the design of synthetic oligonucleotides to mutagenize the non-structural proteins. 

B2.NGS_analysis: Analysis pipeline to get VirVarseq tables, including the modified VirVarSeq script to detect deletions. 

B3.codon_tables_to_MFE: R analysis pipeline to compute MFE from codon tables. 

B4.correction_MFE: R script and data for the correction of MFE based on overlaps between regions and experimental data (competition assay). 

B5.differential_section: R analysis pipeline to compute differential selection between cell lines from codon tables. 