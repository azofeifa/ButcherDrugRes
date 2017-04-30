'''visualize output of DEseq
'''

import pandas as pd,matplotlib.pyplot as plt
import seaborn as sns,numpy as np
def viz(df):
	df 		= df.dropna()
	cp 		= sns.color_palette()
	sns.set(font_scale=2);sns.set_style("ticks")
	colors 	= [cp[2] if row.padj < pow(10,-3) else cp[0]  for i,row in df.iterrows() ]

	F 			= plt.figure(tight_layout=True,figsize=(15,7))
	ax1,ax2 	= F.add_subplot(1,2,1),F.add_subplot(1,2,2)
	ax1.scatter(df["log2FoldChange"],-np.log10(df.pvalue) ,color=colors)
	ax1.set_xlabel("Log 2 Fold Change (A2/AC)");ax1.set_ylabel("-Log 10 P-value")
	ax1.set_title("Volcanco Plot");ax2.set_title("MA Plot")

	ax2.scatter(df["baseMean"],df["log2FoldChange"],color=colors)
	ax2.set_xscale("log")
	ax2.set_ylabel("Log 2 Fold Change (A2/AC)");ax2.set_xlabel("Mean of Expression A2,AC")
	sns.despine()

	plt.show()



def main():
	table 	= "../DEseq/out/Enhancers.csv"
	df 		= pd.read_csv(table)
	viz(df)


main()