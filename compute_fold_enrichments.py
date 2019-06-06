import read_diff_expr_data as data
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import hypergeom
from statsmodels.sandbox.stats.multicomp import multipletests
import math, xlsxwriter

term_order = ['Mitochondrion-localized',
              'Redox-related', 
              ' ']
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

experiments = []
terms = []
fold_enrs = []
ps = []
adjps = []
def compHyperGemP(hits, bg_all, bg_hits, query_size):
    '''
    Returns P(X>=occ);fold-enrichment
    '''
    p = 1-hypergeom.cdf(hits-1, bg_all, bg_hits, query_size) # calculates P(X>=occ)
    fe = hits / ((bg_hits / bg_all) * query_size)
    return p, fe

def checkTerms(): # uses that variables are global
    lexp = spl[0]
    ldir = spl[1]
    
    # compute enrichment
    lps = []
    lfe = []
    for term in term_order:
        if term == " ":
            p, fe = 5.0, 0.0
        else:
            p, fe = compHyperGemP(d_anno[term], bg_s, bg_anno[term], len(d))
        lps.append(p)
        lfe.append(fe)
        
    # adjustment, take care of the fake annotation term
    t_i = lps.index(5.0)
    lps.remove(5.0)
    _, ladjps, _, _ = multipletests(lps, alpha = 0.05, method = "fdr_bh", is_sorted = False, returnsorted = False)
    ladjps = list(ladjps)
    lps.insert(t_i, 1.0)
    ladjps.insert(t_i, 1.0)
    
    # add remaining
    if ldir == "(downregulated)":
        ps.extend(lps)
        fold_enrs.extend(lfe)
        adjps.extend(ladjps)
        experiments.extend( [lexp] * len(term_order))
        terms.extend(term_order)
    
    # build string
    s = [(x, y, z) for x,y,z in zip(lfe, lps, ladjps)]
    return s

data.anno_data[" "] = set() # trick to get an empty line
bg_anno = {}
enr_list = xlsxwriter.Workbook("../excel_out/fold-enr.xlsx")
cell_format = enr_list.add_format({'bg_color': '#33cc33'})

def addListToWorksheet(ws_name, data):
    '''
    Adds list of genes to worksheet and annotates with hallmarks
    '''
    ws = enr_list.add_worksheet(ws_name)
    ws.freeze_panes(1, 0)
    #ws.set_column('B:M', 20)
    ws.write(0, 0, "Annotation term")
    ws.write(0, 1, "fold-enrichment")
    ws.write(0, 2, "p-value")
    ws.write(0, 3, "adj. p-value")
    ws.set_column('A:A', 25)
    ws.set_column('B:B', 15)
    ws.set_column('D:D', 10)
    j = 1
    for term, d in zip(term_order, data):
        fe, pv, apv = d
        form = None
        if term == " ":
            continue
        if apv <= 0.05:
            form = cell_format
        ws.write(j, 0, term, form)
        ws.write_number(j, 1, fe, form)
        ws.write_number(j, 2, pv, form)
        ws.write_number(j, 3, apv, form)
        j += 1
                
for d, s in [(data.nmeth_bg, "Aibar (background)"), (data.nmeth_up, "Aibar (upregulated)"), (data.nmeth_do, "Aibar (downregulated)"),
             (data.gse_bg, "Shoshan (background)"), (data.gse_up, "Shoshan (upregulated)"), (data.gse_do, "Shoshan (downregulated)"), 
             (data.both_bg, "both (background)"), (data.both_up, "both (upregulated)"), (data.both_do, "both (downregulated)")]:
    d_anno = {}
    for term in term_order:
        d_anno[term] = len(d.intersection(data.anno_data[term]))
    
    spl = s.split(" ")
    s2 = spl[0] + "_" + spl[1][1:-1]
    
    # background distr is set
    if "background" in s:
        bg_s = len(d)
        for term in data.anno_data:
            bg_anno[term] = d_anno[term]
        continue
    
    sl = checkTerms()
    addListToWorksheet(s, sl)

enr_list.close()
#
# plot terms of interest
#

sns.set(style="whitegrid")
sns.set_context("talk")
plt.figure(figsize=(10, 7))

# reorder data
idata = pd.DataFrame.from_dict({"experiment":experiments, "annotation term":terms, "fold-enrichment":fold_enrs, "pvalue":ps, "adj. pvalue":adjps})
experiments = []
terms = []
fold_enrs = []
ps = []
adjps = []
logp = []
for term in term_order:
    t = idata[idata["annotation term"] == term]
    for exp in ["Aibar", "Shoshan", "both"]:
        e=t[t["experiment"] == exp]
        experiments.append(exp)
        terms.append(term)
        fold_enrs.append(float(e["fold-enrichment"]))
        ps.append(float(e["pvalue"]))
        p = float(e["adj. pvalue"])
        adjps.append(p)
        if p > 0:
            lp = -math.log10(float(e["adj. pvalue"]))
        else:
            lp = 0
        logp.append(lp)
        
data = pd.DataFrame.from_dict({"experiment":experiments, "annotation term":terms, "fold-enrichment":fold_enrs, "pvalue":ps, "adj. pvalue":adjps, "-log10(adj. p)":logp})

plt.axvline(1, -0.05, 11, color = "crimson", alpha = 0.3)
ax = sns.barplot(y = "annotation term", x = "fold-enrichment", hue = "experiment", data = data, orient = "h", palette = "muted", order = term_order)
plt.legend(frameon = True)
plt.xticks([0, 1, 2, 4, 6, 8, 10])
plt.xlim(0, 9)


def checkP(p):
    s = ""
    if p < 0.05:
        s = "*"
    if p < 0.01:
        s = "**"
    if p < 0.005:
        s = "***"
    return s

def annotateBars(row, ax=ax):
    if row["annotation term"] == " ":
        return
    
    color = "black"
    vertpad = 0.07
    horfactor = 0.334
    horpad = -0.12
    
    upper = range(0, len(data)+1, 3)
    lower = range(2, len(data)+1, 3)
    varpad = 0.0
    varam = 0.06
    if row.name in upper:
        varpad = varam
    elif row.name in lower:
        varpad = -varam
    

    ax.text(row["fold-enrichment"] + vertpad, row.name * horfactor + horpad + varpad, checkP(row["adj. pvalue"]), zorder = 10, rotation = 0, color = color, fontsize = 12, weight = "heavy")
    #ax.text(row["fold-enrichment"] + vertpad, row.name * horfactor + horpad + varpad, row["experiment"] + row["annotation term"], zorder = 10, rotation = 0, color = color, fontsize = 8, weight = "heavy")
    
    return

    vertpad2 = 0.03
    horpad2 = -0.13
    if row["fold-enrichment"] > 0.0:
        ax.text(vertpad2, row.name * horfactor + horpad + varpad + horpad2, "{:1.1f}".format(row["fold-enrichment"]), zorder = 10, rotation = 0, color = color, fontsize = 8, weight = "light")

# invert data for simpler matching in graph
data.reindex(index=data.index[::-1])
data.apply(annotateBars, ax = ax, axis = 1)

plt.tight_layout()
plt.savefig("other_plots/term_fold_enrichment.pdf")