Pythonscripte:
- simulateTrees.py -> used for the simulation of all tree data with the corresponding parameter (parameter were saved in Tree_data.csv)
- graphFunction.py -> all function used for the analysis of the trees and LDT graph (Cluster deletion, get_ldt_triple usw.)
- analyze_plot_tree.py -> Script that iterates over all simulated trees and analyses each tree. The results are appended and saved to Tree_data.csv.

R-Script:
- tree_analyze -> Input: Tree_data.csv. Graphical and non-graphical Evaluation of the data obtained by analyze_plot_tree.py. The output are all plots and csv files.
- graphentheory_Bericht: R-Markdown file to generate the 'Graphentheory_Bericht.pdf' from the Tree_data_Full.csv file.

Plots (in subfolder 02_Plots):
- all plots that were used in the 'Graphentheory_Bericht.pdf'. In the pdf file the plots are described more closly.

Csv-files:
- numeric overview above the used plots
- Sum_Tripple_Fraction -> Tripple_S_Fractions_of_S_Tripples_Groups and Tripple_T_Fractions_of_T_Tripples_Groups
- Result_Gene_vs_HGT -> P0_Gene_vs_HGT to P6_Gene_vs_HGT
- Result_Spezies_vs_HGT -> P0_Spezies_vs_HGT to P6_Spezies_vs_HGT
- Summary_CD_Accuracy -> Accuracy_CD_Groups
- Summary_CD_Precision -> Prec_CD_Groups
- Summary_CD_Recall -> Recall_CD_Groups
- Summary_RS_Accuracy -> Accuracy_RS_Groups
- Summary_RS_Precision -> Prec_RS_Groups
- Summary_RS_Recall -> Recall_RS_Groups
- Tree_data -> all simulation data 
- Tree_data_Full -> all simulation data plus all data used and generated for the evaluation