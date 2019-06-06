from statsmodels.stats.multitest import multipletests
import xlsxwriter

anno_data_file = "../diff_expr_data/annotation_data_iea.tsv"

nmeth_up = set()
nmeth_do = set()
nmeth_bg = set()
nmeth_g2log2fc = {}
for line in open("../diff_expr_data/nmeth_paper.csv"):
    if line.startswith("Gene"):
        continue
    line = line.strip().split(";")
    
    gene = line[0]
    nmeth_bg.add(gene)
    
    # not detected / not expressed
    if line[1] == "0":
        break # -> continue for better background? polyA -> better like this?
    
    log2fc = float(line[2].replace(",", "."))
    nmeth_g2log2fc[gene] = log2fc
             
    if ("_" in line[-1]):
        if "Down" in line[-1]:
            nmeth_do.add(gene)
        elif "UP" in line[-1]:
            nmeth_up.add(gene)

print("dataset", "#genes up", "#genes down", "#genes bg", sep = ";")
print("aibar", len(nmeth_up), len(nmeth_do), len(nmeth_bg), sep = ";")

gse_up = set()
gse_do = set()
gse_bg = set()
gse_p = {}
gse_g2log2fc = {}

#pval adjustment
pvals = []
for line in open("../diff_expr_data/GSE76541_normalized.csv"):
    if line.startswith(";") or line.startswith("ID"):
        continue
    line = line.strip().split(";")
    # ID;KD;NT;KD;NT;KD;NT;KD/NT;KD/NT;Gene_symbol;Description;Group;Entrez_gene;RefSeq;Ensembl_gene;Probe_priority;Probe_number;RNA_nucleotide_gi;Protein_accession.version;Ensembl_transcript;Ensembl_protein
    gene = line[9].split("|")[0] # partially several divided by |
    gse_bg.add(gene)
    
    log2fc = line[7] # log2 ratio KD/NT
    pval = line[8] # pval
    
    if log2fc == "NA" or pval == "NA":
        continue
    
    pval = float(line[8].replace(",", ".")) # pval
    pvals.append(pval)

_, adj_pvals, _, _ = multipletests(pvals, 0.05, "fdr_bh", is_sorted = False, returnsorted = False)

pmap = {}
for p, pad in zip(pvals, adj_pvals):
    pmap[p] = pad
 
# read data   
for line in open("../diff_expr_data/GSE76541_normalized.csv"):
    if line.startswith(";") or line.startswith("ID"):
        continue
    line = line.strip().split(";")
    # ID;KD;NT;KD;NT;KD;NT;KD/NT;KD/NT;Gene_symbol;Description;Group;Entrez_gene;RefSeq;Ensembl_gene;Probe_priority;Probe_number;RNA_nucleotide_gi;Protein_accession.version;Ensembl_transcript;Ensembl_protein
    gene = line[9].split("|")[0] # partially several divided by |
    gse_bg.add(gene)
    
    log2fc = line[7] # log2 ratio KD/NT
    pval = line[8] # pval
    
    if log2fc == "NA" or pval == "NA":
        continue
    
    log2fc = float(line[7].replace(",", ".")) # log2 ratio KD/NT
    gse_g2log2fc[gene] = log2fc
    pval = float(line[8].replace(",", ".")) # pval
    gse_p[gene] = pval
    
    if pmap[pval] < 0.05:# and (log2fc >= 0.2667639314 or log2fc <= -0.2809717645):
        if log2fc > 0:
            gse_up.add(gene)
        elif log2fc < 0:
            gse_do.add(gene)
            
print("shoshan", len(gse_up), len(gse_do), len(gse_bg), sep = ";")

both_up = gse_up.intersection(nmeth_up)
both_do = gse_do.intersection(nmeth_do)
both_bg = gse_bg.intersection(nmeth_bg)

print ("both", len(both_up), len(both_do), len(both_bg), sep = ";")


#
# GO stuff reading
#

def readGOList(p, iea = False):
    genes = set()
    for line in open(p):
        if line.startswith("Symbol"):
            continue
        line = line.strip().split("\t")
        gene, evidence = line
        if gene == "-":
            continue
        gene = gene.upper()
        if evidence == "IEA" and iea == False:
            continue
        genes.add(gene)
        
    return genes

anno_data = {}
for line in open(anno_data_file):
    term, genes = line.strip().split("\t")
    if term == "Calcium ion transmembrane transporter activity":
        continue
    term = term.title()
    if term == "Mitochondrion":
        term = "Mitochondrion-localized"
    if term == "Oxidoreductase Activity" or term == "Antioxidant Activity":
        term = "Redox-related"
    anno_data[term] = set(genes.split(","))

print()
for term in anno_data:
    print(term, len(anno_data[term]), sep = ";")
    
    
#
# write up/down Excel sheet
#

term_order = ['Mitochondrion-localized',
              'Redox-related']
hallmark_terms = ['Deregulating Cellular Energetics', 
                  'Resisting Cell Death', 
                  'Tumor-Promoting Inflammation',
                  'Genome Instability And Mutation',
                  'Evading Growth Suppressors', 
                  'Avoiding Immune Destruction',
                  'Sustaining Proliferative Signaling', 
                  'Activating Invasion And Metastasis', 
              'Enabling Replicative Immortality', 
              'Inducing Angiogenesis']


term_order.extend(hallmark_terms)

up_down_list = xlsxwriter.Workbook("../excel_out/diff_expr_genes.xlsx")
cell_format = up_down_list.add_format({'bg_color': '#33cc33'})

def addListToWorksheet(ws_name, genes):
    '''
    Adds list of genes to worksheet and annotates with hallmarks
    '''
    ws = up_down_list.add_worksheet(ws_name)
    l = sorted([x for x in genes])
    ws.freeze_panes(1, 0)
    ws.set_column('B:M', 20)
    ws.write(0, 0, "Gene")
    for j, term in enumerate(term_order):
        ws.write(0, j+1, term)
    for i, g in enumerate(l):
        ws.write(i+1, 0, g)
        for j, term in enumerate(term_order):
            if g in anno_data[term]:
                ws.write(i+1, j+1, " ", cell_format)
                ws.write(i+1, 0, g, cell_format)
    
for d, s in [(nmeth_up, "Aibar (upregulated)"), (nmeth_do, "Aibar (downregulated)"),  (nmeth_bg, "Aibar (background)"),
             (gse_up, "Shoshan (upregulated)"), (gse_do, "Shoshan (downregulated)"), (gse_bg, "Shoshan (background)"), 
             (both_up, "both (upregulated)"), (both_do, "both (downregulated)"), (both_bg, "both (background)")]:
    
    addListToWorksheet(s, d)

up_down_list.close()