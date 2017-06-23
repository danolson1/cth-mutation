# cth-mutation
Analysis of mutations in engineered strains of Clostridium thermocellum

This repository has a set of tools for analyzing mutations from resequencing data
Here is the basic outline for how to use it:
1. Process raw sequence data with CLC Genomics Workbench to generate CSV files with mutations
2. Run the "process_clc_files_v5.py" script to compile the mutations into a Pandas Dataframe for further analysis.  
3. If desired, the mutations can be annotated with nearby protein coding regions (CDS) by running "annotate_df_with_CDS_v5.py"
4. Finally, a table of mutations can be created by running the "output_pivot_table()" function in "annotate_df_with_CDS_v5.py"
