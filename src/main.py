'''
	changes in motif displacement following drug resistent 
	cell lines treated and not treated with Cisplaton
'''

import sys
sys.path.append("../motif_displacement_analysis/MDAP")
from md_stats import mds_frame
import os,pandas as pd,seaborn as sns,matplotlib.pyplot as plt
from display import show


def main():
	RECOMP= True
	DIR 	= "/Users/joazofeifa/Lab/ButcherDrugRes/mds_files/"

	OUT 	= "../misc_files/A2N_ACN_diff_mds.csv"
	if RECOMP:
		MDS  	= mds_frame()
		f1 	= DIR+"A2N_MDS.csv" 
		f2 	= DIR+"CAN_MDS.csv"
		MDS.load_MD_score_file(f1, "A2N")
		MDS.load_MD_score_file(f2, "ACN")
		df 		= MDS.differential( "A2N", "ACN",h=150)
		df.to_csv(OUT,index=False)

def get_figs():
	DIR 	= "/Users/joazofeifa/Lab/ButcherDrugRes/misc_files/"
	OUT 	= "../images/"
	for f in os.listdir(DIR):
		if "diff_mds" in f:
			print f
			df 			= pd.read_csv(DIR+f)
			sns.set(font_scale=2)
			sns.set_style("ticks")

			F 				= plt.figure(tight_layout=True,figsize=(12,10))	
			ax1,ax2 		= F.add_subplot(1,2,1),F.add_subplot(1,2,2)
			ax2.set_title(f.split(".")[0])
			S 				= show(df=df,ax=ax2,FDR=pow(10,-1.5))
			ax1.hist(df.pval_mds,bins=60,color="steelblue",edgecolor="white",label="N="+str(df.shape[0]))
			ax1.set_xlabel("p-value");ax1.set_ylabel("frequency")
			ax1.legend(loc="best")
			sns.despine()
			plt.savefig(OUT+f.split(".")[0])


main()
get_figs()

