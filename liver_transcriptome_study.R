# Expected expression trends in healthy vs diseased liver
# 150 healthy + 150 diseased (MASLD/MASH, fibrosis F0-F4) liver RNA-seq profiles,
# simulated from a negative-binomial model with biologically structured disease effects,
# then analysed for (1) expression distributions, (2) expected trends / differential
# expression, (3) prediction of disease state and of whole transcriptomes.
#
# Run:  Rscript liver_transcriptome_study.R      (outputs go to ./results)

suppressPackageStartupMessages({
  library(MASS); library(limma); library(DESeq2); library(matrixStats)
})
set.seed(20261006)
dir.create("results", showWarnings = FALSE)
t0 <- Sys.time()
say <- function(...) cat(sprintf(...), "\n", sep = "")

## 1. Cohort ------------------------------------------------------------------
stage_levels <- c("Healthy", "F0-F1", "F2", "F3", "F4")
stage <- factor(c(rep("Healthy", 150), rep("F0-F1", 50), rep("F2", 40),
                  rep("F3", 35), rep("F4", 25)), levels = stage_levels)
n <- length(stage)
k <- as.integer(stage) - 1                      # 0 = healthy ... 4 = F4
group <- factor(ifelse(k == 0, "Healthy", "Diseased"), levels = c("Healthy", "Diseased"))
H <- group == "Healthy"; D <- !H
sev <- ifelse(k == 0, 0, pmin(pmax((k + rnorm(n, 0, 0.35)) / 4, 0.05), 1.15))  # latent severity
sex <- factor(sample(c("F", "M"), n, replace = TRUE))
batch <- factor(sample(c("B1", "B2"), n, replace = TRUE))
lib_target <- round(20e6 * exp(rnorm(n, 0, 0.3)))

## 2. Genes and modules -------------------------------------------------------
# amp = log2 fold change at full severity (s = 1); shape = trend over severity s in [0, 1.15]
mods <- list(
  metabolic = list(label = "Hepatocyte metabolism (CYP)", n = 200, amp = -1.4, shape = function(s) s,
    named = c(CYP2E1 = 12, CYP3A4 = 10.5, CYP2C9 = 10.5, CYP2C8 = 10, CYP1A2 = 9.5, GNMT = 8.5,
              CYP4A11 = 8, SLC22A1 = 8, G6PC1 = 8, CYP8B1 = 7.5, CYP2C19 = 7, HAO2 = 7), mu = 6, sd = 2),
  lipid = list(label = "Lipid handling / steatosis", n = 80, amp = 1.5, shape = function(s) sin(pi * pmin(s, 1) * 0.8),
    named = c(PLIN2 = 6.5, SCD = 6, FASN = 5, SREBF1 = 5, PNPLA3 = 4.5, ACACA = 4.5, ELOVL6 = 3.5,
              CIDEC = 3, AKR1B10 = 2.5, CIDEA = 0.5), mu = 4, sd = 1.5),
  inflam = list(label = "Inflammation / macrophage", n = 120, amp = 1.6, shape = function(s) s^0.7,
    named = c(CD74 = 8.5, `HLA-DRA` = 7.5, LYZ = 7, CD68 = 6.5, CD9 = 5, GPNMB = 4, CCL2 = 4,
              FABP5 = 4, CXCL10 = 3, TREM2 = 1.5), mu = 4, sd = 1.5),
  ecm = list(label = "Fibrosis / extracellular matrix", n = 150, amp = 2.6, shape = function(s) s^1.5,
    named = c(TIMP1 = 6, COL1A2 = 5.5, COL3A1 = 5.5, COL1A1 = 5, LUM = 5, ACTA2 = 5, EFEMP1 = 4.5,
              MMP2 = 4.5, FBLN5 = 4, PDGFRB = 3.5, LOX = 2.5, THBS2 = 2), mu = 3, sd = 1.5),
  ductular = list(label = "Ductular reaction", n = 60, amp = 2.6, shape = function(s) s^3,
    named = c(SPP1 = 4, CD24 = 3.5, KRT19 = 3, SOX9 = 3, KRT7 = 2.5, EPCAM = 2), mu = 1.5, sd = 1.5),
  synthetic = list(label = "Plasma-protein synthesis", n = 60, amp = -0.9, shape = function(s) s^2,
    named = c(ALB = 15.5, APOA1 = 14, APOC3 = 13, FGB = 12.5, AHSG = 12, HPX = 12, TF = 12,
              RBP4 = 11.5, TTR = 11, ITIH4 = 10.5), mu = 9, sd = 2)
)
G <- 6000
gene <- sprintf("G%04d", seq_len(G)); module <- rep("background", G)
base <- ifelse(runif(G) < 0.3, rnorm(G, -1, 1.5), rnorm(G, 4, 2.2))   # log2 relative abundance
w <- numeric(G); pos <- 0
for (m in names(mods)) {
  md <- mods[[m]]; ids <- pos + seq_len(md$n); pos <- pos + md$n
  module[ids] <- m
  base[ids] <- rnorm(md$n, md$mu, md$sd)
  w[ids] <- ifelse(runif(md$n) < 0.15, runif(md$n, 0.05, 0.3), pmax(rnorm(md$n, 1, 0.25), 0.3))
  nn <- length(md$named)
  gene[ids[1:nn]] <- names(md$named); base[ids[1:nn]] <- md$named
  w[ids[1:nn]] <- pmax(rnorm(nn, 1.05, 0.15), 0.6)
}
sex_genes <- c(XIST = "F", RPS4Y1 = "M", DDX3Y = "M", KDM5D = "M", UTY = "M", EIF1AY = "M")
sids <- pos + seq_along(sex_genes); pos <- pos + length(sex_genes)
gene[sids] <- names(sex_genes); module[sids] <- "sex"; base[sids] <- 4
bg <- which(module == "background")
bg_de <- sample(bg, round(0.05 * length(bg)))            # 5% of background genes get small effects
bg_amp <- numeric(G); bg_amp[bg_de] <- rnorm(length(bg_de), 0, 0.5)

## 3. Expected log2 expression and counts ------------------------------------
lfc <- matrix(0, G, n)
for (m in names(mods)) {
  ids <- which(module == m)
  lfc[ids, ] <- outer(mods[[m]]$amp * w[ids], mods[[m]]$shape(sev))
}
lfc[bg_de, ] <- outer(bg_amp[bg_de], sev)
Fac <- matrix(rnorm(n * 4), n, 4); L <- matrix(rnorm(G * 4, 0, 0.12), G, 4)   # shared latent factors
tau <- ifelse(H, 0.15, 0.35)                                                  # module activity variation
U <- sapply(names(mods), function(m) rnorm(n) * tau)
U[, "ecm"] <- 0.6 * U[, "inflam"] + 0.8 * U[, "ecm"]
modrand <- matrix(0, G, n)
for (m in names(mods)) { ids <- which(module == m); modrand[ids, ] <- outer(w[ids], U[, m]) }
batch_eff <- outer(rnorm(G, 0, 0.1), as.numeric(batch == "B2"))
sex_eff <- matrix(0, G, n)
for (j in seq_along(sids)) sex_eff[sids[j], ] <- ifelse(sex == sex_genes[j], 0, -8)
eta <- base + lfc + L %*% t(Fac) + modrand + batch_eff + sex_eff
p <- 2^eta; p <- t(t(p) / colSums(p))
mu <- t(t(p) * lib_target)
phi <- exp(rnorm(G, log(0.03), 0.6))
phi_mat <- outer(phi, ifelse(H, 1, 1.2))
counts <- matrix(rnbinom(G * n, mu = mu, size = 1 / phi_mat), G, n,
                 dimnames = list(gene, sprintf("S%03d", seq_len(n))))
true_lfc <- rowMeans(lfc[, D])
nonnull <- module %in% names(mods) | seq_len(G) %in% bg_de
meta <- data.frame(sample = colnames(counts), group, stage, severity = round(sev, 3), sex, batch,
                   library_size = colSums(counts))
write.csv(meta, "results/sample_metadata.csv", row.names = FALSE)
gz <- gzfile("results/liver_counts.csv.gz", "w"); write.csv(counts, gz); close(gz)
write.csv(data.frame(gene, module, base_log2 = round(base, 3), true_mean_log2FC = round(true_lfc, 4),
                     nb_dispersion = round(phi, 4)), "results/gene_truth.csv", row.names = FALSE)

## 4. Expression distributions ------------------------------------------------
sink("results/summary.txt", split = TRUE)
say("== Cohort ==")
print(table(stage, sex))
say("Library size: median %.1f M (range %.1f-%.1f M)", median(meta$library_size) / 1e6,
    min(meta$library_size) / 1e6, max(meta$library_size) / 1e6)

cd <- data.frame(group, sex, batch, row.names = colnames(counts))
dds <- DESeqDataSetFromMatrix(counts, cd, ~ batch + sex + group)
dds <- estimateSizeFactors(dds); sf <- sizeFactors(dds)
lib <- colSums(counts)
ncpm <- t(t(counts) / (sf * mean(lib))) * 1e6
logcpm <- log2(ncpm + 0.5)
nc <- t(t(counts) / sf)

say("\n== Distribution ==")
for (g in c("Healthy", "Diseased")) {
  s <- group == g
  det <- colSums(counts[, s] >= 5)
  zero <- colMeans(counts[, s] == 0)
  top10 <- apply(counts[, s], 2, function(x) sum(sort(x, decreasing = TRUE)[1:10]) / sum(x))
  alb <- counts["ALB", s] / lib[s]
  gm <- rowMeans(logcpm[, s])
  say("%s: genes >=5 counts median %d; zero fraction %.3f; top10 share %.3f; ALB share %.3f; gene mean log2CPM quantiles 10/25/50/75/90 = %s",
      g, as.integer(median(det)), median(zero), median(top10), median(alb),
      paste(round(quantile(gm, c(.1, .25, .5, .75, .9)), 2), collapse = "/"))
}
expressed <- rowMeans(logcpm) > 1
say("Expressed genes (mean log2CPM > 1): %d of %d", sum(expressed), G)

disp_group <- function(s) {
  m <- rowMeans(nc[, s]); v <- rowVars(nc[, s])
  phi_hat <- (v - m * mean(1 / sf[s])) / m^2
  sqrt(pmax(phi_hat, 0))
}
bcv_H <- disp_group(H); bcv_D <- disp_group(D)
amean <- rowMeans(logcpm)
bins <- cut(amean, c(-Inf, 0, 2, 4, 6, 8, 10, Inf), labels = c("<0", "0-2", "2-4", "4-6", "6-8", "8-10", ">10"))
disp_bins <- do.call(rbind, lapply(levels(bins), function(b) {
  ii <- bins == b
  data.frame(bin = b, genes = sum(ii), bcv_healthy = round(median(bcv_H[ii]), 3),
             bcv_diseased = round(median(bcv_D[ii]), 3))
}))
print(disp_bins); write.csv(disp_bins, "results/dispersion_by_expression.csv", row.names = FALSE)
say("True NB BCV (healthy) median: %.3f; ratio diseased/healthy BCV over expressed genes: %.2f",
    median(sqrt(phi)), median(bcv_D[expressed] / bcv_H[expressed]))

vH <- rowVars(logcpm[, H]); vD <- rowVars(logcpm[, D])
Fst <- vD / vH
pF <- 2 * pmin(pf(Fst, 149, 149), 1 - pf(Fst, 149, 149))
qF <- p.adjust(pF[expressed], "BH")
say("Variance test (expressed genes): more variable in disease %d, less variable %d (FDR<0.05)",
    sum(qF < 0.05 & Fst[expressed] > 1), sum(qF < 0.05 & Fst[expressed] < 1))
# share of diseased-group variance explained by fibrosis stage
xD <- logcpm[, D]; stD <- droplevels(stage[D])
stage_means <- sapply(levels(stD), function(s) rowMeans(xD[, stD == s]))
resid <- xD - stage_means[, as.integer(stD)]
within <- rowSums(resid^2) / (sum(D) - nlevels(stD))
explained <- 1 - within / vD
modset <- c(names(mods), "background")
mod_tab <- do.call(rbind, lapply(modset, function(m) {
  ii <- module == m & expressed
  data.frame(module = m, genes = sum(ii), bcv_healthy = round(median(bcv_H[ii]), 3),
             bcv_diseased = round(median(bcv_D[ii]), 3),
             var_ratio_D_over_H = round(median(Fst[ii]), 2),
             share_more_variable = round(mean(qF[ii[expressed]] < 0.05 & Fst[expressed][ii[expressed]] > 1), 3),
             stage_explained = round(median(explained[ii]), 3))
}))
print(mod_tab); write.csv(mod_tab, "results/module_variance.csv", row.names = FALSE)

## 5. Differential expression and expected trends ----------------------------
say("\n== Differential expression (DESeq2, ~ batch + sex + group) ==")
dds <- DESeq(dds, quiet = TRUE)
res <- results(dds, contrast = c("group", "Diseased", "Healthy"))
tested <- !is.na(res$padj)
sig <- tested & res$padj < 0.05
say("Tested %d genes; significant %d (up %d, down %d); |log2FC|>1 & padj<0.05: %d (up %d, down %d)",
    sum(tested), sum(sig), sum(sig & res$log2FoldChange > 0), sum(sig & res$log2FoldChange < 0),
    sum(sig & abs(res$log2FoldChange) > 1), sum(sig & res$log2FoldChange > 1), sum(sig & res$log2FoldChange < -1))
say("Empirical FDR %.3f; power for |true log2FC|>=0.5: %.3f; power for |true log2FC| 0.1-0.5: %.3f",
    sum(sig & !nonnull) / sum(sig),
    sum(sig & nonnull & abs(true_lfc) >= 0.5) / sum(tested & nonnull & abs(true_lfc) >= 0.5),
    sum(sig & nonnull & abs(true_lfc) >= 0.1 & abs(true_lfc) < 0.5) /
      sum(tested & nonnull & abs(true_lfc) >= 0.1 & abs(true_lfc) < 0.5))
ok <- tested & nonnull & expressed
say("Correlation estimated vs true log2FC (non-null expressed genes, n=%d): r = %.3f; slope = %.3f",
    sum(ok), cor(res$log2FoldChange[ok], true_lfc[ok]), coef(lm(res$log2FoldChange[ok] ~ true_lfc[ok]))[2])
for (m in names(mods)) {
  ii <- module == m & tested
  say("  %-10s DE %3d/%3d (%.0f%%), median est log2FC %.2f (true %.2f)", m, sum(sig[ii]), sum(ii),
      100 * mean(sig[ii]), median(res$log2FoldChange[ii]), median(true_lfc[ii]))
}
markers <- c("COL1A1", "COL3A1", "THBS2", "LUM", "TREM2", "GPNMB", "CXCL10", "AKR1B10", "PLIN2",
             "KRT7", "SPP1", "CYP2C19", "CYP1A2", "CYP2E1", "ALB", "APOC3")
mk <- data.frame(gene = markers, module = module[match(markers, gene)],
                 healthy_log2CPM = round(rowMeans(logcpm[markers, H]), 1),
                 true_log2FC = round(true_lfc[match(markers, gene)], 2),
                 est_log2FC = round(res[markers, "log2FoldChange"], 2),
                 padj = signif(res[markers, "padj"], 2))
print(mk); write.csv(mk, "results/marker_genes.csv", row.names = FALSE)
res_out <- data.frame(gene, module, baseMean = round(res$baseMean, 1), log2FC = round(res$log2FoldChange, 3),
                      lfcSE = round(res$lfcSE, 3), pvalue = signif(res$pvalue, 3), padj = signif(res$padj, 3),
                      true_log2FC = round(true_lfc, 3))
write.csv(res_out[order(res_out$padj), ], "results/deseq2_diseased_vs_healthy.csv", row.names = FALSE)

say("\n== Module trends across stages (module score = mean log2 change vs healthy mean) ==")
hmean <- rowMeans(logcpm[, H])
trend <- do.call(rbind, lapply(names(mods), function(m) {
  ii <- module == m & expressed
  score <- colMeans(logcpm[ii, ] - hmean[ii])
  expct <- colMeans(lfc[ii, ])
  do.call(rbind, lapply(stage_levels, function(s) {
    sc <- score[stage == s]
    data.frame(module = m, label = mods[[m]]$label, stage = s, median = round(median(sc), 3),
               q1 = round(quantile(sc, .25), 3), q3 = round(quantile(sc, .75), 3),
               expected = round(median(expct[stage == s]), 3))
  }))
}))
rownames(trend) <- NULL
print(trend); write.csv(trend, "results/module_stage_trend.csv", row.names = FALSE)
say("Observed vs expected stage medians: r = %.3f, mean abs diff = %.3f",
    cor(trend$median, trend$expected), mean(abs(trend$median - trend$expected)))

say("\n== limma trend test: log2CPM ~ orthogonal poly(stage, 2) + sex + batch ==")
P2 <- poly(k, 2)
des <- model.matrix(~ P2 + sex + batch)
fit <- eBayes(lmFit(logcpm[expressed, ], des), trend = TRUE)
q_lin <- p.adjust(fit$p.value[, "P21"], "BH"); q_quad <- p.adjust(fit$p.value[, "P22"], "BH")
b_lin <- fit$coef[, "P21"]; b_quad <- fit$coef[, "P22"]
say("Expressed genes with linear trend %d, with curvature %d (FDR<0.05)", sum(q_lin < 0.05), sum(q_quad < 0.05))
mod_e <- module[expressed]
curv <- do.call(rbind, lapply(names(mods), function(m) {
  ii <- mod_e == m
  data.frame(module = m, genes = sum(ii),
             pct_linear_up = round(100 * mean(q_lin[ii] < 0.05 & b_lin[ii] > 0)),
             pct_linear_down = round(100 * mean(q_lin[ii] < 0.05 & b_lin[ii] < 0)),
             pct_convex = round(100 * mean(q_quad[ii] < 0.05 & b_quad[ii] > 0)),
             pct_concave = round(100 * mean(q_quad[ii] < 0.05 & b_quad[ii] < 0)))
}))
print(curv); write.csv(curv, "results/module_curvature.csv", row.names = FALSE)
bgq <- mod_e == "background"
say("Background genes with curvature: %.1f%%", 100 * mean(q_quad[bgq] < 0.05))

## 6. Prediction --------------------------------------------------------------
X <- t(logcpm[expressed, ]); y <- D
make_folds <- function(seed) {
  set.seed(seed); f <- integer(n)
  for (s in stage_levels) { ii <- which(stage == s); f[ii] <- sample(rep(1:5, length.out = length(ii))) }
  f
}
auc <- function(score, y) { r <- rank(score); n1 <- sum(y); n0 <- sum(!y); (sum(r[y]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
ridge_core <- function(Xtr, Ytr) {
  mx <- colMeans(Xtr); sx <- colSds(Xtr); sx[sx == 0] <- 1
  Z <- sweep(sweep(Xtr, 2, mx), 2, sx, "/")
  my <- colMeans(Ytr); s <- svd(Z)
  list(mx = mx, sx = sx, my = my, s = s, Uty = crossprod(s$u, sweep(Ytr, 2, my)))
}
# Ridge regression; one penalty shared by all outputs, chosen by inner 5-fold CV
# (GCV breaks down when predictors outnumber samples, so it is not used).
ridge_fit <- function(Xtr, Ytr, lambdas = 10^seq(-2, 6, length.out = 33), inner = 5) {
  Ytr <- as.matrix(Ytr); set.seed(99); fi <- sample(rep(1:inner, length.out = nrow(Xtr)))
  sse <- numeric(length(lambdas))
  for (f in 1:inner) {
    a <- fi != f; ca <- ridge_core(Xtr[a, , drop = FALSE], Ytr[a, , drop = FALSE])
    ZV <- sweep(sweep(Xtr[!a, , drop = FALSE], 2, ca$mx), 2, ca$sx, "/") %*% ca$s$v
    Yb <- sweep(Ytr[!a, , drop = FALSE], 2, ca$my)
    for (j in seq_along(lambdas))
      sse[j] <- sse[j] + sum((Yb - ZV %*% ((ca$s$d / (ca$s$d^2 + lambdas[j])) * ca$Uty))^2)
  }
  jl <- which.min(sse); l <- lambdas[jl]
  if (jl %in% c(1, length(lambdas))) say("  note: lambda %.3g at the edge of the grid", l)
  cr <- ridge_core(Xtr, Ytr)
  list(mx = cr$mx, sx = cr$sx, my = cr$my, B = cr$s$v %*% ((cr$s$d / (cr$s$d^2 + l)) * cr$Uty), lambda = l)
}
ridge_predict <- function(fit, Xte) {
  Z <- sweep(sweep(Xte, 2, fit$mx), 2, fit$sx, "/")
  sweep(Z %*% fit$B, 2, fit$my, "+")
}

say("\n== Disease classification, 10 x 5-fold CV ==")
cls <- list(); det_stage <- matrix(0, 10, 5, dimnames = list(NULL, stage_levels))
for (r in 1:10) {
  folds <- make_folds(100 + r)
  p_pca <- p_sig <- p_one <- numeric(n); c_pca <- c_sig <- c_one <- logical(n)
  for (f in 1:5) {
    tr <- folds != f; te <- !tr
    mtr <- colMeans(X[tr, ]); str <- colSds(X[tr, ]); str[str == 0] <- 1
    Xc_tr <- sweep(X[tr, ], 2, mtr); Xc_te <- sweep(X[te, ], 2, mtr)
    v <- svd(Xc_tr, nu = 0, nv = 10)$v
    Ztr <- Xc_tr %*% v; Zte <- Xc_te %*% v
    ld <- lda(Ztr, grouping = y[tr]); pr <- predict(ld, Zte)
    p_pca[te] <- pr$posterior[, "TRUE"]; c_pca[te] <- pr$class == "TRUE"
    # 20-gene signature (10 up, 10 down) chosen on the training fold only
    m1 <- colMeans(X[tr & y, ]); m0 <- colMeans(X[tr & !y, ])
    sp <- sqrt((colVars(X[tr & y, ]) + colVars(X[tr & !y, ])) / 2) + 1e-6
    tt <- (m1 - m0) / sp
    up <- order(tt, decreasing = TRUE)[1:10]; dn <- order(tt)[1:10]
    zsc <- function(idx) (rowMeans(sweep(sweep(X[idx, up, drop = FALSE], 2, mtr[up]), 2, str[up], "/")) -
                          rowMeans(sweep(sweep(X[idx, dn, drop = FALSE], 2, mtr[dn]), 2, str[dn], "/")))
    s_tr <- zsc(which(tr)); s_te <- zsc(which(te))
    g1 <- suppressWarnings(glm(y[tr] ~ s_tr, family = binomial))
    p_sig[te] <- s_te; c_sig[te] <- (coef(g1)[1] + coef(g1)[2] * s_te) > 0
    # single gene: COL1A1
    x1 <- X[, "COL1A1"]
    g2 <- glm(y[tr] ~ x1[tr], family = binomial)
    p_one[te] <- x1[te]; c_one[te] <- (coef(g2)[1] + coef(g2)[2] * x1[te]) > 0
  }
  for (nm in c("pca", "sig", "one")) {
    pp <- get(paste0("p_", nm)); cc <- get(paste0("c_", nm))
    cls[[length(cls) + 1]] <- data.frame(model = nm, rep = r, auc = auc(pp, y), acc = mean(cc == y),
                                         sens = mean(cc[y]), spec = mean(!cc[!y]))
  }
  det_stage[r, ] <- tapply(c_pca, stage, mean)
}
cls <- do.call(rbind, cls)
cls_sum <- aggregate(cbind(auc, acc, sens, spec) ~ model, cls, function(x) round(mean(x), 3))
cls_sd <- aggregate(auc ~ model, cls, function(x) round(sd(x), 3))
print(merge(cls_sum, cls_sd, by = "model", suffixes = c("", "_sd")))
write.csv(cls, "results/classification_cv.csv", row.names = FALSE)
say("Called diseased by stage (PCA-LDA, mean over repeats): %s",
    paste(sprintf("%s %.1f%%", stage_levels, 100 * colMeans(det_stage)), collapse = "; "))

say("\n== Fibrosis-stage regression (ridge on all expressed genes, 5-fold CV) ==")
folds <- make_folds(7); pk <- numeric(n); lam <- c()
for (f in 1:5) {
  tr <- folds != f; te <- !tr
  rf <- ridge_fit(X[tr, ], k[tr]); pk[te] <- ridge_predict(rf, X[te, ]); lam <- c(lam, rf$lambda)
}
say("Spearman rho %.3f; within one stage %.1f%%; exact stage %.1f%%; lambdas %s",
    cor(pk, k, method = "spearman"), 100 * mean(abs(round(pmin(pmax(pk, 0), 4)) - k) <= 1),
    100 * mean(round(pmin(pmax(pk, 0), 4)) == k), paste(signif(lam, 2), collapse = ","))
say("Mean predicted stage by true stage: %s",
    paste(sprintf("%s %.2f", stage_levels, tapply(pk, stage, mean)), collapse = "; "))

say("\n== Transcriptome prediction (5-fold CV) ==")
set.seed(11)
eg <- colnames(X)
eval_genes <- sample(eg, floor(length(eg) / 2)); pool <- setdiff(eg, eval_genes)
say("Evaluation genes %d; landmark pool %d", length(eval_genes), length(pool))
Ks <- c(10, 25, 50, 100, 200, 400)
models <- c("Global mean", "Disease label", "Fibrosis stage", paste0("LM", Ks), paste0("LM", Ks, "+stage"))
SSE <- matrix(0, length(eval_genes), length(models), dimnames = list(eval_genes, models)); SST <- SSE
samp_r <- matrix(NA, n, length(models), dimnames = list(NULL, models))
folds <- make_folds(21)
stage_dm <- model.matrix(~ stage)[, -1]
for (f in 1:5) {
  tr <- folds != f; te <- !tr
  Ytr <- X[tr, eval_genes]; Yte <- X[te, eval_genes]; mtr <- colMeans(Ytr)
  preds <- list()
  preds[["Global mean"]] <- matrix(mtr, sum(te), length(eval_genes), byrow = TRUE)
  gm_ <- rbind(Healthy = colMeans(Ytr[y[tr] == FALSE, ]), Diseased = colMeans(Ytr[y[tr] == TRUE, ]))
  preds[["Disease label"]] <- gm_[ifelse(y[te], "Diseased", "Healthy"), ]
  sm_ <- t(sapply(stage_levels, function(s) colMeans(Ytr[stage[tr] == s, ])))
  preds[["Fibrosis stage"]] <- sm_[as.character(stage[te]), ]
  Zp <- scale(X[tr, pool])
  for (K in Ks) {
    set.seed(1000 + K + f)
    km <- suppressWarnings(kmeans(t(Zp), centers = K, iter.max = 50))
    lm_genes <- sapply(seq_len(K), function(cl) {       # gene closest to each cluster centre
      ii <- which(km$cluster == cl)
      ii[which.min(colSums((Zp[, ii, drop = FALSE] - km$centers[cl, ])^2))]
    })
    lmg <- pool[lm_genes]
    rf <- ridge_fit(X[tr, lmg], Ytr); preds[[paste0("LM", K)]] <- ridge_predict(rf, X[te, lmg])
    rf2 <- ridge_fit(cbind(X[tr, lmg], stage_dm[tr, ]), Ytr)
    preds[[paste0("LM", K, "+stage")]] <- ridge_predict(rf2, cbind(X[te, lmg], stage_dm[te, ]))
    say("  fold %d: K=%d lambda %.3g, +stage lambda %.3g", f, K, rf$lambda, rf2$lambda)
  }
  for (md in models) {
    P <- preds[[md]]
    SSE[, md] <- SSE[, md] + colSums((Yte - P)^2)
    SST[, md] <- SST[, md] + colSums(sweep(Yte, 2, mtr)^2)
    dev_obs <- sweep(Yte, 2, mtr); dev_pred <- sweep(P, 2, mtr)
    samp_r[which(te), md] <- sapply(seq_len(sum(te)), function(i)
      if (sd(dev_pred[i, ]) == 0) 0 else cor(dev_obs[i, ], dev_pred[i, ]))
  }
}
R2 <- 1 - SSE / SST
pred_sum <- data.frame(model = models, median_gene_R2 = round(apply(R2, 2, median), 3),
                       pct_genes_R2_gt_0.5 = round(100 * colMeans(R2 > 0.5), 1),
                       median_sample_r = round(apply(samp_r, 2, median), 3),
                       median_sample_r_healthy = round(apply(samp_r[H, ], 2, median), 3),
                       median_sample_r_diseased = round(apply(samp_r[D, ], 2, median), 3))
print(pred_sum, row.names = FALSE); write.csv(pred_sum, "results/transcriptome_prediction.csv", row.names = FALSE)
# oracle ceiling: squared correlation of observed log2CPM with the true expected log2 abundance,
# shifted by each sample's normalisation offset (median over expressed genes)
off <- apply(logcpm[expressed, ] - eta[expressed, ], 2, median)
oracle <- sapply(eval_genes, function(g) cor(X[, g], eta[match(g, gene), ] + off)^2)
say("Oracle R2 (true expected expression), median over evaluation genes: %.3f", median(oracle))
emod <- module[match(eval_genes, gene)]
bymod <- do.call(rbind, lapply(c(names(mods), "background"), function(m) {
  ii <- emod == m
  data.frame(module = m, genes = sum(ii), R2_label = round(median(R2[ii, "Disease label"]), 3),
             R2_stage = round(median(R2[ii, "Fibrosis stage"]), 3),
             R2_LM100 = round(median(R2[ii, "LM100"]), 3), R2_LM100_stage = round(median(R2[ii, "LM100+stage"]), 3),
             R2_oracle = round(median(oracle[ii]), 3))
}))
print(bymod); write.csv(bymod, "results/transcriptome_prediction_by_module.csv", row.names = FALSE)
# trivial (uncentred) profile correlation, for contrast
say("Uncentred profile correlation, observed vs global-mean prediction (median over samples): %.4f",
    median(sapply(seq_len(n), function(i) cor(X[i, eval_genes], colMeans(X[, eval_genes])))))
say("\nRuntime %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins")))
sink()
