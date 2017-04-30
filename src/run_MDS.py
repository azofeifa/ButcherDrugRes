'''
	changes in motif displacement following drug resistent 
	cell lines treated and not treated with Cisplaton
'''

import sys,numpy as np
sys.path.append("../motif_displacement_analysis/MDAP")
sys.path.append("../pyD3/")

from md_stats import mds_frame
import os,pandas as pd,seaborn as sns,matplotlib.pyplot as plt
from display import show
import HTML
F 		= "A2C_ACC"

def calc():
	RECOMP= True
	DIR 	= "../mds_files/ntss/"
	OUT 	= "../mds_files/ordered_pvalues/" +F+ "_diff_mds.csv"

	MDS  	= mds_frame()
	f1 	= DIR+F.split("_")[0]+"_MDS.csv" 
	f2 	= DIR+F.split("_")[1]+"_MDS.csv"
	MDS.load_MD_score_file(f1, F.split("_")[0])
	MDS.load_MD_score_file(f2, F.split("_")[1])
	df 		= MDS.differential( F.split("_")[0], F.split("_")[1],h=150)
	df.to_csv(OUT,index=False)

def viz():
	df 	= pd.read_csv("../mds_files/ordered_pvalues/" + F + "_diff_mds.csv").dropna()
	ax 	= HTML.axes()
	ax.scatter(list(df.mds),-np.log10(df.pval_mds+pow(10,-15)),alpha=1.0,size=3,lbls=[x.split("_")[1]for x in list(df.motif)])
	ax.set_xlabel("MDS Difference (" + F.split("_")[1] + "-" + F.split("_")[0] + ")")
	ax.set_ylabel("-log 10 p-value")
	ax.savefig("../images/" + F + ".html")
	ax.show()

calc()
viz()