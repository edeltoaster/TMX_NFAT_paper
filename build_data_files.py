import os, gzip, ClinicalData
import analyze_BRAF_V600E_dep as braf

cd = ClinicalData.ClinicalData("../clinical/clinical.tsv")

mutect2 = braf.getBRAFstate("../mutation_data/TCGA.SKCM.mutect2.somatic.maf.gz")

def readExpr(p):
    '''
    reads expression file
    '''
    data = {}
    for line in gzip.open(p):
        line = line.decode("utf-8").strip().split()
        gene = line[0]
        expr = float(line[1])
        data[gene] = expr
    return data

surv_data = cd.getSurvivalData() # patient -> (time in days, status)
stages = cd.getAnnotation("tumor_stage")
ages = cd.getAge()
genders = cd.getGender()
to_write = []
to_write.append(["patient", "age", "sex", "status", "survival (days)", "survival (years)", "NFAT1 expr", "TMX1 expr", "TMX3 expr", "tumor stage", "BRAF"])
to_write_R = []
to_write_R.append(["patient", "age", "sex", "status", "survival_days", "survival_years", "NFAT1_expr", "TMX1_expr", "TMX3_expr", "tumor_stage", "BRAF"])
for fi in os.listdir("../SKCM/"):
    if not fi.endswith(".txt.gz"):
        continue
    p = fi.split(".")[0]
    time, status = surv_data[p]
    time_y = float(time) / 365.25
    expr = readExpr("../SKCM/" + fi)
    nfat1_expr = expr["NFATC2"]
    tmx1_expr = expr["TMX1"]
    tmx3_expr = expr["TMX3"]
    stage = stages[p]
    age = ages[p]
    gender = genders[p]
    # "patient", "age", "sex", "status", "survival (days)", "survival (years)", "NFAT1 expr", "TMX1 expr", "TMX3 expr", "tumor stage"
    to_write.append([p, str(age), gender, status, str(time), str(time_y), str(nfat1_expr), str(tmx1_expr), str(tmx3_expr), stage, mutect2.get(p, "WT")])
    if status == "True":
        status = "TRUE"
    else:
        status = "FALSE"
    to_write_R.append([p, str(age), gender, status, str(time), str(time_y), str(nfat1_expr), str(tmx1_expr), str(tmx3_expr), stage, mutect2.get(p, "WT")])

open("relevant_SKCM_data.csv", "w").writelines([";".join(x)+"\n" for x in to_write])
open("relevant_SKCM_data_R.csv", "w").writelines([";".join(x)+"\n" for x in to_write_R])