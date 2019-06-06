#install.packages("ggplot2")
#install.packages("survminer")
#install.packages("ggfortify")
library(survival)
library(survminer)
library(stringr)

setwd("~/GDrive/Work/projects/NFATC2_TMX/scripts")
df <- read.csv("relevant_SKCM_data_R.csv", na.strings = "NA", sep = ";")
mask <- df["survival_days"]>0 # clean for negative surv patient
df <- df[mask,]
model <- Surv(df$survival_years, df$status)
df["status_col"] <- ifelse(df["status"] == TRUE, "red", "blue")
df["tumor_stage_num"] <- ifelse(df[,"tumor_stage"] == "stage iv", 4, str_count(df[,"tumor_stage"], "i"))

#
# descriptive stuff
#

summary(df)

cor(df$NFAT1_expr,df$TMX1_expr) # 0.17, not correlated
cor(df$NFAT1_expr,df$TMX3_expr) # 0.18 not correlated
cor(df$TMX1_expr,df$TMX3_expr) # 0.67

# not correlated
cor(df$NFAT1_expr, df$tumor_stage_num) # 0.1
cor(df$TMX1_expr, df$tumor_stage_num) # -0.05
cor(df$TMX3_expr, df$tumor_stage_num) # 0.12

#
# optimize thresholds
#

to_test_nfat <- sort(c(df[,"NFAT1_expr"]))
q_nfat <- quantile(to_test_nfat, c(0.2, 0.8))
mask_nfat <- (q_nfat[1] < to_test_nfat) & (to_test_nfat < q_nfat[2])
to_test_nfat <- to_test_nfat[mask_nfat]

to_test_tmx1 <- sort(c(df[,"TMX1_expr"]))
q_tmx1 <- quantile(to_test_tmx1, c(0.2, 0.8))
mask_tmx1 <- (q_tmx1[1] < to_test_tmx1) & (to_test_tmx1 < q_tmx1[2])
to_test_tmx1 <- to_test_tmx1[mask_tmx1]

to_test_tmx3 <- sort(c(df[,"TMX3_expr"]))
q_tmx3 <- quantile(to_test_tmx3, c(0.2, 0.8))
mask_tmx3 <- (q_tmx3[1] < to_test_tmx3) & (to_test_tmx3 < q_tmx3[2])
to_test_tmx3 <- to_test_tmx3[mask_tmx3]

to_test_age <- sort(c(df[,"age"]))
q_age <- quantile(to_test_age, c(0.2, 0.8))
mask_age <- (q_age[1] < to_test_age) & (to_test_age < q_age[2])
to_test_age <- to_test_age[mask_age]

nfat_opt_p <- 1
nfat_opt_thr <- 1
for (th_nfat in to_test_nfat) {
  df["NFAT1_state"] <- c()
  df["NFAT1_state"] <- ifelse(df["NFAT1_expr"] > th_nfat, "high", "low")
  mht <- survdiff(model~df$NFAT1_state)
  pval <- 1 - pchisq(mht$chisq, length(mht$n) - 1)
  
  if (pval < nfat_opt_p) {
    nfat_opt_p <- pval
    nfat_opt_thr <- th_nfat
  }
}
df["NFAT1_state"] <- ifelse(df["NFAT1_expr"] > nfat_opt_thr, 1, 0)
print(paste0("NFAT1 thr:", nfat_opt_thr, " p:", nfat_opt_p))

tmx1_opt_p <- 1
tmx1_opt_thr <- 1
for (th_tmx1 in to_test_tmx1) {
  df["TMX1_state"] <- c()
  df["TMX1_state"] <- ifelse(df["TMX1_expr"] > th_tmx1, "high", "low")
  mht <- survdiff(model~df$TMX1_state)
  pval <- 1 - pchisq(mht$chisq, length(mht$n) - 1)
  
  if (pval < tmx1_opt_p) {
    tmx1_opt_p <- pval
    tmx1_opt_thr <- th_tmx1
  }
}
df["TMX1_state"] <- ifelse(df["TMX1_expr"] > tmx1_opt_thr, 1, 0)
print(paste0("TMX1 thr:", tmx1_opt_thr, " p:", tmx1_opt_p))

tmx3_opt_p <- 1
tmx3_opt_thr <- 1
for (th_tmx3 in to_test_tmx3) {
  df["TMX3_state"] <- c()
  df["TMX3_state"] <- ifelse(df["TMX3_expr"] > th_tmx3, "high", "low")
  mht <- survdiff(model~df$TMX3_state)
  pval <- 1 - pchisq(mht$chisq, length(mht$n) - 1)
  
  if (pval < tmx3_opt_p) {
    tmx3_opt_p <- pval
    tmx3_opt_thr <- th_tmx3
  }
}
df["TMX3_state"] <- ifelse(df["TMX3_expr"] > tmx3_opt_thr, 1, 0)
print(paste0("TMX3 thr:", tmx3_opt_thr, " p:", tmx3_opt_p))

age_opt_p <- 1
age_opt_thr <- 1
for (th_age in to_test_age) {
  df["age_state"] <- c()
  df["age_state"] <- ifelse(df["age"] > th_age, "older", "younger")
  mht <- survdiff(model~df$age_state)
  pval <- 1 - pchisq(mht$chisq, length(mht$n) - 1)
  
  if (pval < age_opt_p) {
    age_opt_p <- pval
    age_opt_thr <- th_age
  }
}
df["age_state"] <- ifelse(df["age"] > age_opt_thr, 1, 0)
print(paste0("age thr:", age_opt_thr, " p:", age_opt_p))

pdf("other_plots/TMX1_stage_boxplot.pdf")
boxplot(df$TMX1_expr ~ df$tumor_stage, las = 2, ylab = "TMX1 expression [FPKM]")
abline(h = tmx1_opt_thr, col = "tomato") 
dev.off()
pdf("other_plots/TMX3_stage_boxplot.pdf")
boxplot(df$TMX3_expr ~ df$tumor_stage, las = 2, ylab = "TMX3 expression [FPKM]")
abline(h = tmx3_opt_thr, col = "tomato") 
dev.off()
pdf("other_plots/NFAT1_stage_boxplot.pdf")
boxplot(df$NFAT1_expr ~ df$tumor_stage, las = 2, ylab = "NFAT1 expression [FPKM]")
abline(h = nfat_opt_thr, col = "tomato")
dev.off()

#
# state stubs
#
df[df$NFAT1_state = 0 & df$BRAF]
df["NFAT1_age_state"] <- paste(df[,"NFAT1_state"], df[,"age_state"], sep = "/")
df["TMX1_age_state"] <- paste(df[,"TMX1_state"], df[,"age_state"], sep = "/")
df["TMX3_age_state"] <- paste(df[,"TMX3_state"], df[,"age_state"], sep = "/")
df["NFAT1_sex_state"] <- paste(df[,"NFAT1_state"], df[,"sex"], sep = "/")
df["TMX1_sex_state"] <- paste(df[,"TMX1_state"], df[,"sex"], sep = "/")
df["TMX3_sex_state"] <- paste(df[,"TMX3_state"], df[,"sex"], sep = "/")

df["exact_state"] <- paste(df[,"NFAT1_state"], df[,"TMX1_state"], df[,"TMX3_state"], sep = "/")
df["simple_state"] <- (df[,"NFAT1_state"] == 1) + (df[,"TMX1_state"] == 1) + (df[,"TMX3_state"]== 1)
df["bool_mintwo"] <- ifelse( (df[,"simple_state"] > 1), 1, 0)
df["bool_any"] <- ifelse( (df[,"simple_state"] > 0), 1, 0)
df["NT3_pair"] <- ifelse((df[,"NFAT1_state"] == 1) & (df[,"TMX3_state"] == 1), 1, 0)
df["NT1_pair"] <- ifelse((df[,"NFAT1_state"] == 1) & (df[,"TMX1_state"] == 1), 1, 0)
df["T13_pair"] <- ifelse((df[,"TMX1_state"] == 1) & (df[,"TMX3_state"] == 1), 1, 0)
df["triplet"] <- ifelse((df[,"NFAT1_state"] == 1) & (df[,"TMX1_state"] == 1) & (df[,"TMX3_state"] == 1), 1, 0)
df["simple_BRAF"] <- ifelse(df[,"BRAF"] == "WT", 0, 1)

nf_thrs <- quantile(df[,"NFAT1_expr"], c(0.33, 0.66))
df["NFAT1_tert"] <- ifelse( df[,"NFAT1_expr"] < nf_thrs[1], 0, ifelse( df[,"NFAT1_expr"] < nf_thrs[2], 1, 2))
tmx1_thrs <- quantile(df[,"TMX1_expr"], c(0.33, 0.66))
df["TMX1_tert"] <- ifelse( df[,"TMX1_expr"] < tmx1_thrs[1], 0, ifelse( df[,"TMX1_expr"] < tmx1_thrs[2], 1, 2))
tmx3_thrs <- quantile(df[,"TMX3_expr"], c(0.33, 0.66))
df["TMX3_tert"] <- ifelse( df[,"TMX3_expr"] < tmx3_thrs[1], 0, ifelse( df[,"TMX3_expr"] < tmx3_thrs[2], 1, 2))

#
# until here to load general data
#


#
# conditional probabiities ?
#
ggplot(df, aes(x=TMX3_expr, fill=sex), ggtheme = theme_minimal()) + geom_density(alpha=.2)
ggplot(df, aes(x=TMX1_expr, fill=sex), ggtheme = theme_minimal()) + geom_density(alpha=.2)
ggplot(df, aes(x=NFAT1_expr, fill=sex), ggtheme = theme_minimal()) + geom_density(alpha=.2)


#
# Kaplan Meier and log-rank test
#

descriptors <- colnames(df)
mask <- c("patient", "age", "status", "survival_days", "survival_years", "NFAT1_expr", "TMX1_expr", "TMX3_expr", "status_col", "tumor_stage")
descriptors <- descriptors [! descriptors %in% mask]

# automatic checking
descrs <- c()
Ps <- c()
for (groups in descriptors) {
  model <- Surv(df$survival_years, df$status)
  f <- as.formula(paste0("model ~ df$", groups))
  mht <- do.call(survdiff, args = list(formula = f, data = df))
  pval <- 1 - pchisq(mht$chisq, length(mht$n) - 1)
  print(paste0(groups, " p:", pval))
  descrs <- append(descrs, groups)
  Ps <- append(Ps, pval)
  kme <- do.call(survfit, args = list(formula = f, data = df))
  plots <- ggsurvplot(kme, data = df, xlab = "Time [years]", risk.table = TRUE, pval = TRUE,conf.int = TRUE, conf.int.alpha = 0.15, ggtheme = theme_minimal())
  ggsave(filename = paste0("surv_plots_all/", groups, ".pdf"), print(plots, newpage = FALSE), height = 7, width = 7)
}
lr_results <- data.frame(descrs, Ps)
colnames(lr_results) <- c("descriptor", "log-rank pval")
write.csv(lr_results, "surv_outs/descr_all.csv", row.names = FALSE, quote = FALSE)

# subset check, more than 1 worthwhile
mask <- df["simple_state"]>0 # any increase
df2 <- df[mask,]
descrs <- c()
Ps <- c()
for (groups in descriptors) {
  if (length(unique(df2[,groups])) == 1) {
    next
  }
  model <- Surv(df2$survival_years, df2$status)
  f <- as.formula(paste0("model ~ df2$", groups))
  mht <- do.call(survdiff, args = list(formula = f, data = df2)) # could also be done using subset = ..., same result anyhow
  pval <- 1 - pchisq(mht$chisq, length(mht$n) - 1)
  print(paste0(groups, " p:", pval))
  descrs <- append(descrs, groups)
  Ps <- append(Ps, pval)
  kme <- do.call(survfit, args = list(formula = f, data = df2))
  plots <- ggsurvplot(kme, data = df2, xlab = "Time [years]", risk.table = TRUE, pval = TRUE,conf.int = TRUE, conf.int.alpha = 0.15, ggtheme = theme_minimal())
  ggsave(filename = paste0("surv_plots_any/", groups, ".pdf"), print(plots, newpage = FALSE), height = 7, width = 7)
}
lr_results <- data.frame(descrs, Ps)
colnames(lr_results) <- c("descriptor", "log-rank pval")
write.csv(lr_results, "surv_outs/descr_any.csv", row.names = FALSE, quote = FALSE)

# check of sex/age/synergy effect on NFAT1, TMX1 and TMX3
descrs <- c()
Ps <- c()
for (goi in c("NFAT1_state", "TMX1_state", "TMX3_state")) {
  mask <- df[goi] > 0 # gene of interest increased
  df2 <- df[mask,]
  for (groups in c("age_state", "sex", "NFAT1_state", "TMX1_state", "TMX3_state")) {
    if (length(unique(df2[,groups])) == 1) {
      next
    }
    model <- Surv(df2$survival_years, df2$status)
    f <- as.formula(paste0("model ~ df2$", groups))
    mht <- do.call(survdiff, args = list(formula = f, data = df2)) # could also be done using subset = ..., same result anyhow
    pval <- 1 - pchisq(mht$chisq, length(mht$n) - 1)
    print(paste0(goi, "-", groups, " p:", pval))
    descrs <- append(descrs, paste0(goi, "-", groups))
    Ps <- append(Ps, pval)
    kme <- do.call(survfit, args = list(formula = f, data = df2))
    plots <- ggsurvplot(kme, data = df2, xlab = "Time [years]", risk.table = TRUE, pval = TRUE,conf.int = TRUE, conf.int.alpha = 0.15, ggtheme = theme_minimal())
    ggsave(filename = paste0("surv_plots_specific/", goi, "-", groups, ".pdf"), print(plots, newpage = FALSE), height = 7, width = 7)
  }  
}
lr_results <- data.frame(descrs, Ps)
colnames(lr_results) <- c("descriptor", "log-rank pval")
write.csv(lr_results, "surv_outs/descr_specific.csv", row.names = FALSE, quote = FALSE)


#
# Cox Prop. Hazard models
#

df["mult_expr"] <- df[, "NFAT1_expr"] * df[, "TMX1_expr"] * df[, "TMX3_expr"]
df["NT1_mult"] <- df[, "NFAT1_expr"] * df[, "TMX1_expr"]
df["NT3_mult"] <- df[, "NFAT1_expr"] * df[, "TMX3_expr"]
df["T13_mult"] <- df[, "TMX1_expr"] * df[, "TMX3_expr"]

descriptors <- colnames(df)
mask <- c("patient", "status", "survival_days", "survival_years", "status_col", "tumor_stage")
descriptors <- descriptors [! descriptors %in% mask]

# univariate
descrs <- c()
HRs <- c()
CI_ls <- c()
CI_hs <- c()
zscores <- c()
Prs <- c()
lr_tests <- c()
model <- Surv(df$survival_years, df$status)
for (descr in descriptors) {
  f <- as.formula(paste0("model ~ df$", descr))
  cox <- do.call(coxph, args = list(formula = f, data = df))
  scox <- summary(cox)
  HR <- scox$coefficients[2]
  CI_l <- scox$conf.int[3]
  CI_h <- scox$conf.int[4]
  zscore <- scox$coefficients[4]
  Pr <- scox$coefficients[5]
  lr_test_p <- scox$logtest[3]
  print(paste(descr, HR, CI_l, CI_h, zscore, Pr, lr_test_p, sep = " , "))
  descrs <- append(descrs, descr)
  HRs <- append(HRs, HR)
  CI_ls <- append(CI_ls, CI_l)
  CI_hs <- append(CI_hs, CI_h)
  zscores <- append(zscores, zscore)
  Prs <- append(Prs, Pr)
  lr_tests <- append(lr_tests, lr_test_p)
}
univariate_cox_results <- data.frame(descrs, HRs, CI_ls, CI_hs, zscores, Prs, lr_tests)
colnames(univariate_cox_results) <- c("descriptor", "HR", "CI low", "CI high", "z-Score", "pval", "LR-test pval")
write.csv(univariate_cox_results, "surv_outs/univariate_cox_all.csv", row.names = FALSE, quote = FALSE)


#
# multivariate COX
#

model <- Surv(df$survival_years, df$status)
summary(coxph(model ~ df$NFAT1_expr + df$TMX1_expr + df$TMX3_state + df$age + df$sex, data = df))
summary(coxph(model ~ df$NFAT1_expr + df$TMX1_expr + df$TMX3_state + df$age + df$sex, data = df))
summary(coxph(model ~ df$NFAT1_expr + df$TMX3_state + df$age + df$sex, data = df))
summary(coxph(model ~ df$NFAT1_expr + df$TMX3_state, data = df))
summary(coxph(model ~ df$simple_state + df$age + df$sex, data = df))
summary(coxph(model ~ df$simple_state, data = df))
summary(coxph(model ~ df$NFAT1_expr + df$TMX3_state, data = df))
summary(coxph(model ~ df$NFAT1_expr + df$TMX1_expr + df$TMX3_state + df$age_state + df$sex, data = df))
summary(coxph(model ~ df$NFAT1_expr + df$TMX3_state+ df$NT3_mult, data = df))
summary(coxph(model ~ df$NFAT1_state + df$TMX3_state, data = df))
summary(coxph(model ~ df$TMX3_state, data = df))

# nice explanation: http://www.sthda.com/english/wiki/cox-model-assumptions#testing-proportional-hazards-assumption

cox <- (coxph(model ~ df$NFAT1_expr + df$TMX3_state, data = df))
testcox <- cox.zph(cox) # test for proportional hazard assumption -> ok if nothing significant
plot <- ggcoxzph(testcox) # scaled Schoenfeld residuals against the transformed time tested for each covariate
ggsave(filename = paste0("other_plots/schoenfeld_residues.pdf"), print(plot, newpage = FALSE), height = 7, width = 7)

ggcoxdiagnostics(cox, type = "schoenfeld", ox.scale = "time")
ggcoxdiagnostics(cox, type = "dfbeta", linear.predictions = FALSE, ggtheme = theme_bw())
ggcoxdiagnostics(cox, type = "deviance", linear.predictions = FALSE, ggtheme = theme_bw())


#
# for paper
#

palette <- c('#4878CF', '#D65F5F') #c("#E7B800", "#2E9FDF")
ggtheme <- theme_minimal()
alpha <- 0.2
theight <- 0.2

model <- Surv(df$survival_years, df$status)

leg <- c("other patients", "NFAT1+")
kme <- survfit(model ~ df$NFAT1_state) 
plots <- ggsurvplot(kme, data = df, xlab = "Time [years]", risk.table = TRUE, pval = TRUE,conf.int = TRUE, conf.int.alpha = alpha, ggtheme = ggtheme, palette = palette, tables.height = theight, legend.labs = leg)
ggsave(filename = paste0("other_plots/NFAT1.pdf"), print(plots, newpage = FALSE), height = 7, width = 7)

leg <- c("other patients", "TMX1+")
kme <- survfit(model ~ df$TMX1_state) 
plots <- ggsurvplot(kme, data = df, xlab = "Time [years]", risk.table = TRUE, pval = TRUE,conf.int = TRUE, conf.int.alpha = alpha, ggtheme = ggtheme, palette = palette, tables.height = theight, legend.labs = leg)
ggsave(filename = paste0("other_plots/TMX1.pdf"), print(plots, newpage = FALSE), height = 7, width = 7)

leg <- c("other patients", "TMX3+")
kme <- survfit(model ~ df$TMX3_state) 
plots <- ggsurvplot(kme, data = df, xlab = "Time [years]", risk.table = TRUE, pval = TRUE,conf.int = TRUE, conf.int.alpha = alpha, ggtheme = ggtheme, palette = palette, tables.height = theight, legend.labs = leg)
ggsave(filename = paste0("other_plots/TMX3.pdf"), print(plots, newpage = FALSE), height = 7, width = 7)

leg <- c("other patients", "any GOI+")
kme <- survfit(model ~ df$bool_any) 
plots <- ggsurvplot(kme, data = df, xlab = "Time [years]", risk.table = TRUE, pval = TRUE,conf.int = TRUE, conf.int.alpha = alpha, ggtheme = ggtheme, palette = palette, tables.height = theight, legend.labs = leg)
ggsave(filename = paste0("other_plots/any_GOI.pdf"), print(plots, newpage = FALSE), height = 7, width = 7)

# univariate Cox
descriptors <- c("NFAT1_state", "TMX1_state", "TMX3_state", "NFAT1_expr", "TMX1_expr", "TMX3_expr")
descrs <- c()
HRs <- c()
CI_ls <- c()
CI_hs <- c()
zscores <- c()
Prs <- c()
lr_tests <- c()
model <- Surv(df$survival_years, df$status)
for (descr in descriptors) {
  f <- as.formula(paste0("model ~ df$", descr))
  cox <- do.call(coxph, args = list(formula = f, data = df))
  scox <- summary(cox)
  HR <- scox$coefficients[2]
  CI_l <- scox$conf.int[3]
  CI_h <- scox$conf.int[4]
  zscore <- scox$coefficients[4]
  Pr <- scox$coefficients[5]
  lr_test_p <- scox$logtest[3]
  print(paste(descr, HR, CI_l, CI_h, zscore, Pr, lr_test_p, sep = " , "))
  descrs <- append(descrs, descr)
  HRs <- append(HRs, HR)
  CI_ls <- append(CI_ls, CI_l)
  CI_hs <- append(CI_hs, CI_h)
  zscores <- append(zscores, zscore)
  Prs <- append(Prs, Pr)
  lr_tests <- append(lr_tests, lr_test_p)
}
univariate_cox_results <- data.frame(descrs, HRs, CI_ls, CI_hs, zscores, Prs, lr_tests)
colnames(univariate_cox_results) <- c("descriptor", "HR", "CI 95%-low", "CI 95%-high", "z-score", "pval", "LR-test pval")
write.csv(univariate_cox_results, "surv_outs/univariate_cox.csv", row.names = FALSE, quote = FALSE)

# -> best non-correlated
cox <- (coxph(model ~ df$NFAT1_expr + df$TMX3_state, data = df))
summary(cox)
testcox <- cox.zph(cox) # test for proportional hazard assumption -> ok if nothing significant
plot <- ggcoxzph(testcox, ggtheme = ggtheme, palette = palette) # scaled Schoenfeld residuals against the transformed time tested for each covariate
ggsave(filename = paste0("other_plots/schoenfeld_residues.pdf"), print(plot, newpage = FALSE), height = 7, width = 7)

ggcoxdiagnostics(cox, type = "schoenfeld", ox.scale = "time")
ggcoxdiagnostics(cox, type = "dfbeta", linear.predictions = FALSE, ggtheme = ggtheme, palette = palette)
ggcoxdiagnostics(cox, type = "deviance", linear.predictions = FALSE, ggtheme = ggtheme, palette = palette)


# for revision
wt <- df[which(df$BRAF == "WT"), ]
mut <- df[which(df$BRAF == "V600E"), ]

prop.table(table(wt$NFAT1_state))*nrow(wt) # 34 low, 15 high
prop.table(table(mut$NFAT1_state))*nrow(mut) # 22 low, 26 high

m <- matrix(c(34,22,15,26),2,2)
dimnames(m) <-  list(c("WT", "V600E"), c("low", "high"))
fisher.test(m) # p=0.02421 , fisher's exact test the better test here

prop.table(table(wt$TMX1_state))*nrow(wt) # 38 low, 11 high
prop.table(table(mut$TMX1_state))*nrow(mut) # 28 low, 20 high

m <- matrix(c(38,28,11,20),2,2)
dimnames(m) <-  list(c("WT", "V600E"), c("low", "high"))
fisher.test(m) # p=0.05169

prop.table(table(wt$TMX3_state))*nrow(wt) # 30 low, 19 high
prop.table(table(mut$TMX3_state))*nrow(mut) # 25 low, 23 high

m <- matrix(c(30,25,19,23),2,2)
dimnames(m) <-  list(c("WT", "V600E"), c("low", "high"))
fisher.test(m) # p=0.42

# check for BRAF (V600E) dependence
descrs <- c()
Ps <- c()
for (state in c("WT", "V600E")) {
  mask <- df["BRAF"] == state # gene of interest increased
  df2 <- df[mask,]
  for (groups in c("NFAT1_state", "TMX1_state", "TMX3_state")) {
    if (length(unique(df2[,groups])) == 1) {
      next
    }
    model <- Surv(df2$survival_years, df2$status)
    f <- as.formula(paste0("model ~ df2$", groups))
    mht <- do.call(survdiff, args = list(formula = f, data = df2)) # could also be done using subset = ..., same result anyhow
    pval <- 1 - pchisq(mht$chisq, length(mht$n) - 1)
    print(paste0(state, "-", groups, " p:", pval))
    descrs <- append(descrs, paste0(state, "-", groups))
    Ps <- append(Ps, pval)
    kme <- do.call(survfit, args = list(formula = f, data = df2))
    plots <- ggsurvplot(kme, data = df2, xlab = "Time [years]", risk.table = TRUE, pval = TRUE,conf.int = TRUE, conf.int.alpha = 0.15, ggtheme = theme_minimal())
    ggsave(filename = paste0("surv_plots_BRAF/", state, "-", groups, ".pdf"), print(plots, newpage = FALSE), height = 7, width = 7)
  }  
}
lr_results <- data.frame(descrs, Ps)
colnames(lr_results) <- c("descriptor", "log-rank pval")
write.csv(lr_results, "surv_outs/BRAF_specific.csv", row.names = FALSE, quote = FALSE)

cox <- (coxph(model ~ df$NFAT1_expr + df$TMX3_state + df$simple_BRAF, data = df))
testcox <- cox.zph(cox) # test for proportional hazard assumption -> ok if nothing significant

#
# materials for revision
#

# new survival plots
palette <- c('#4878CF', '#D65F5F') #c("#E7B800", "#2E9FDF")
ggtheme <- theme_minimal()
alpha <- 0.2
theight <- 0.2

mask <- df["BRAF"] == "WT"
df_wt <- df[mask,]
model <- Surv(df_wt$survival_years, df_wt$status)
leg <- c("other patients", "NFAT1+")
kme <- survfit(model ~ df_wt$NFAT1_state) 
plots <- ggsurvplot(kme, data = df_wt, xlab = "Time [years]", risk.table = TRUE, pval = TRUE,conf.int = TRUE, conf.int.alpha = alpha, ggtheme = ggtheme, palette = palette, tables.height = theight, legend.labs = leg)
ggsave(filename = paste0("other_plots/NFAT1_in_WT.pdf"), print(plots, newpage = FALSE), height = 7, width = 7)

mask <- df["BRAF"] == "V600E"
df_m <- df[mask,]
model <- Surv(df_m$survival_years, df_m$status)
leg <- c("other patients", "NFAT1+")
kme <- survfit(model ~ df_m$NFAT1_state) 
plots <- ggsurvplot(kme, data = df_m, xlab = "Time [years]", risk.table = TRUE, pval = TRUE,conf.int = TRUE, conf.int.alpha = alpha, ggtheme = ggtheme, palette = palette, tables.height = theight, legend.labs = leg)
ggsave(filename = paste0("other_plots/NFAT1_in_V600E.pdf"), print(plots, newpage = FALSE), height = 7, width = 7)

