'''
	take in as input one of the count .tsv tables;
	visualize as heatmap and cluster
'''
import pandas as pd,matplotlib.pyplot as plt
import numpy as np,scipy.stats as ss
import seaborn as sns
import matplotlib as mpl,os
import scipy.cluster.hierarchy as sch
import math
def viz(df,K=500,title="RefSeq Gene Annotation"):
	try:
		df 		= df[["chrom","start","stop","gene", 
									"A2N", "A2C", "A2D", "ACN", "ACC", "ACD"]]
		X 			= np.array(df.as_matrix()[:,4:],dtype=float)
		xlabels 	= df.columns[3:]
	except:
		df 		= df[["chrom","start","stop", 
									"A2N", "A2C", "A2D", "ACN", "ACC", "ACD"]]
		X 			= np.array(df.as_matrix()[:,3:],dtype=float)
		xlabels 	= df.columns[2:]


	lengths 	= df.stop-df.start
	'''filter out non zeros
	'''
	X 			= X[(X > 1).all(axis=1),:]
	'''normalize by length
	'''
	X 			= np.array([X[i,:]/float(lengths[i]) for i in range(X.shape[0])])
	'''log 10 scale
	'''
	X 			= np.log10(X)
	'''Z-score
	'''
	Z 			= ss.zscore(X,axis=0, ddof=1)
	stds 		= np.argsort(np.std(Z,axis=1)/np.mean(Z,axis=1))[::-1][:K]
	stds 		= np.argsort(np.std(Z,axis=1))[::-1][:K]
	Z 			= Z[stds,:]
	'''cluster etc.
	'''
	Y 			= sch.linkage(Z, method='centroid')
	Z1 		= sch.dendrogram(Y, orientation='right',no_plot=True)
	Z 			= Z[Z1["leaves"],:]

	'''draw
	'''
	sns.set(font_scale=1.4); sns.set_style("ticks")
	cmap1 = mpl.colors.ListedColormap(sns.diverging_palette(240, 10))
	F 			= plt.figure(figsize=(5,10),tight_layout=True)
	ax 		= plt.gca()
	ax.set_title(title)
	plt.imshow(Z,vmin=-1.5,vmax=1.5,interpolation="nearest",aspect="auto",cmap=cmap1)
	ax.set_xticklabels(xlabels,rotation=90)
	ax.set_ylabel("top " + str(K)+ " most variable")
	sns.despine()
	cb 	= plt.colorbar(aspect=10,label='Expression (z-score)')
	cb.outline.set_visible(False)
	plt.show()

def load_matrix(G,FILE):
	lines = open(FILE,"r").readlines()[1:]
	for line in lines:
		line_array 	= line.strip("\n").split(",")
		motif,d 		= line_array[0],map(float, line_array[1:])
		if motif not in G:
			G[motif] = list()
		if d[4] > 0:
			EX 			= d[2]*(d[1]/d[3])

			G[motif]+=[d[0]/EX,d[-1]]
		else:
			G[motif]+=[0,d[-1]]
	return G
def viz_ME(G,xid):
	A,pv,M 	= list(),list(),list()
	for motif in G:
		A.append( [G[motif][i] for i in range(len(G[motif])) if i%2==0  ])
		pv.append([G[motif][i] for i in range(len(G[motif])) if i%2==1  ])
		M.append(motif.split("_")[1])
	A 		= np.array(A);pv 	= np.array(pv);M=np.array(M)
	vals 	= np.min(pv,axis=1) < pow(10,-1) 
	vals2 = np.max(A,axis=1) > 5
	A 		= A[vals & vals2]
	M 		= M[vals & vals2][::-1]
	sns.set(font_scale=1.5)
	F 		= plt.figure(figsize=(5,10),tight_layout=True)
	ax 	= plt.gca()
	sns.heatmap(A,edgecolor="white",linewidth=1,ax=ax,cbar_kws={"label":"Fold Change Above Expectation"})
	ax.set_xticklabels([ i.split("_")[1] for i in xid],rotation=90)
	ax.set_yticklabels(M,rotation=0)
	plt.savefig("/Users/joazofeifa/BBC/HDAC8/images/TF_enr_table")
	plt.show()




def main():
	GENES 	= False
	ENHANCE 	= False
	MOTIF_ENR= True
	if MOTIF_ENR:
		G 	= {}
		DIR= "../motif_enrichment/out/"
		xid= list()
		for f in os.listdir(DIR):
			G 	= load_matrix(G,DIR+f)
			xid.append(f)
		viz_ME(G,xid)

	if GENES:
		table 	= "../tables/gene_count_table.tsv"
		title 	= "RefSeq Gene Annotation"		
	if ENHANCE:
		table 	= "../tables/Gm12878_StrongEnhancer_count_table.tsv"
		title 	= "Strong Enhancers\n(Chrom HMM)"
	df 		= pd.read_csv(table,sep="\t")

	viz(df,title=title)

main()