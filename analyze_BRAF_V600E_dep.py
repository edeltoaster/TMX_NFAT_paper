import gzip, os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats 

sns.set(style="whitegrid")
sns.set_context("talk")
plt.figure(figsize=(10, 7))

def getBRAFstate(maf_file):
    
    patient_to_braf = {}
    for line in gzip.open(maf_file):
        line = line.decode("utf8").strip()
        
        if not line.startswith("BRAF"):
            continue
        
        line = line.split("\t")
        
        mutation = line[36]
        braf = ""
        if mutation != "p.V600E":
            braf = "other"
        else:
            braf = "V600E"
        
        patient = "-".join(line[15].split("-")[:3])
        if patient in patient_to_braf and patient_to_braf[patient] == "V600E":
            continue
        patient_to_braf[patient] = braf
        
    return patient_to_braf

def readExprData():
    expr_data = {}
    for line in open("relevant_SKCM_data_R.csv"):
        if line.startswith("patient"):
            continue
        line = line.split(";")
        patient = line[0]
        nfat = float(line[6])
        tmx1 = float(line[7])
        tmx3 = float(line[8])
        expr_data[patient] = (nfat, tmx1, tmx3)
        
    return expr_data

if __name__ == "__main__":
    #muse = getBRAFstate("../mutation_data/TCGA.SKCM.muse.somatic.maf.gz")
    mutect2 = getBRAFstate("../mutation_data/TCGA.SKCM.mutect2.somatic.maf.gz")
    
    expr_data = readExprData()
    patients = []
    mutations = []
    genes = []
    expressions = []
    
    # patient mutation gene expr
    for patient in expr_data:
        patients.extend( [patient]*3 )
        mut = mutect2.get(patient, "WT")
        mutations.extend( [mut]*3 )
        
        genes.append("NFAT1")
        expressions.append(expr_data[patient][0])
        genes.append("TMX1")
        expressions.append(expr_data[patient][1])
        genes.append("TMX3")
        expressions.append(expr_data[patient][2])
        
    d = {"Patient": patients, "BRAF genotype": mutations, "gene": genes, "expression level": expressions}
    
    plot_data = pd.DataFrame(data = d)
    
    print("sizes: WT V600E other")
    wt = plot_data[plot_data["BRAF genotype"] == "WT"]
    v600e = plot_data[plot_data["BRAF genotype"] == "V600E"]
    other = plot_data[plot_data["BRAF genotype"] == "other"]
    print(len(wt)/3, len(v600e)/3, len(other)/3)
    
    print("Mann-Whitney U")
    for g in ["NFAT1", "TMX1", "TMX3"]:
        wt = plot_data[plot_data["BRAF genotype"] == "WT"]
        wt = wt[wt["gene"] == g]["expression level"]
        mutated = plot_data[plot_data["BRAF genotype"] == "V600E"]
        mutated = mutated[mutated["gene"] == g]["expression level"]
        p = stats.mannwhitneyu(wt, mutated, alternative = "two-sided")[1]
        print(g, p)
    
    plot_data = plot_data[plot_data["BRAF genotype"] != "other"]
    sns.boxplot(x = "gene", y = "expression level", hue = "BRAF genotype", data = plot_data, showfliers = False,  palette = "muted")
    plt.tight_layout()
    plt.savefig("other_plots/BRAF_expression_boxplot.pdf")
    plt.clf()

    sns.boxplot(x = "gene", y = "expression level", hue = "BRAF genotype", data = plot_data, showfliers = True,  palette = "muted")
    plt.tight_layout()
    plt.savefig("other_plots/BRAF_expression_boxplot_with_outliers.pdf")
    