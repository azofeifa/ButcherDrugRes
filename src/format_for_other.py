
import time,pandas as pd

def DEseq(FILE,OUT,header,REF=False):
	G 	= dict()
	FHW= open(OUT, "w")
	FHW.write(header)
	i 	= 0
	lines 	= open(FILE,"r").readlines()[1:]
	for line in lines:
		i+=1

		line_array=line.strip("\n").split("\t")
		if line_array[3] not in G and REF:
			G[line_array[3]]=None
			FHW.write("\t".join(line_array[3:])+"\n")
		elif not REF:
			FHW.write(line_array[0]+":"+line_array[1] + "-" + line_array[2] + "\t" + "\t".join(line_array[3:])+"\n")
def overlap(F,OUT):
	df 		= pd.read_csv(F).dropna()
	threshold= pow(10,-3)
	FHW 		= open(OUT, "w")
	for i,row in df.iterrows():
		chrom,stsp 	= i.split(":")
		start,stop 	= stsp.split("-")
		if row.padj < threshold and row["log2FoldChange"] < 0:
			FHW.write(chrom+"\t" + start + "\t" + stop + "\t1\n")
		else:
			FHW.write(chrom+"\t" + start + "\t" + stop + "\t0\n")			



def main():
	REF 	= False
	ENH 	= False
	OV 	= True
	if OV:
		F 	= "../DEseq/out/Enhancers.csv"
		OUT = "../motif_enrichment/files/Sig_Genes_AC.bed"
		overlap(F,OUT)
	if ENH:
		FILE 	= "../tables/Gm12878_StrongEnhancer_count_table.tsv"
		OUT 	= "../DEseq/files/Enhancer_count_table_A2vAC.tsv"
		DEseq(FILE, OUT, "\tA21\tA22\tAC1\tAC2\tAC3\tA23\n", REF=False)		
	if REF:
		FILE 	= "../tables/gene_count_table.tsv"
		OUT 	= "../DEseq/files/RefSeq_count_table_A2vAC.tsv"
		DEseq(FILE, OUT, "\tA21\tA22\tAC1\tAC2\tAC3\tA23\n", REF=True)

main()
