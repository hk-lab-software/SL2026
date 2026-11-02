# =============================================================================
# Extended Fig 4. c
# Hypergeometric over-representation analysis (ORA)
# Cornification / keratinization gene sets + SPRR family
# =============================================================================
setwd("ExtendedFig4c")
## ---- 0. Settings -----------------------------------------------------------
deseq_csv <-"DESeq2_results.csv"   # <-- edit
out_prefix<-"ORA_cornification"
lfc_cut<-1      # |log2FoldChange| threshold
fdr_cut<-0.1    # padj threshold
direction<-"up"   # "up" = log2FC >= lfc_cut; "both" = |log2FC| >= lfc_cut


## ---- 1. Gene sets ----------------------------------------------------------
library(data.table)
GOBP_CORNIFICATION<-fread("Gene_sets/GOBP_CORNIFICATION.v2026.1.Hs.gmt", data.table = F, header = F)
GOBP_CORNIFICATION<-as.character(GOBP_CORNIFICATION[1,-c(1:2)])

GOBP_KERATINIZATION<-fread("Gene_sets/GOBP_KERATINIZATION.v2026.1.Hs.gmt", data.table = F, header = F)
GOBP_KERATINIZATION<-as.character(GOBP_KERATINIZATION[1,-c(1:2)])

REACTOME_CORNIFIED_ENVELOPLE<-fread("Gene_sets/REACTOME_FORMATION_OF_THE_CORNIFIED_ENVELOPE.v2026.1.Hs.gmt",
                                    data.table = F, header = F)
REACTOME_CORNIFIED_ENVELOPLE<-as.character(REACTOME_CORNIFIED_ENVELOPLE[1,-c(1:2)])

EPIDERMAL_DIFFERENTIATION_COMPLEX<-fread("Gene_sets/Epidermal_differentiation_complex.gmt",
                                         data.table = F, header = F)
EPIDERMAL_DIFFERENTIATION_COMPLEX<-as.character(EPIDERMAL_DIFFERENTIATION_COMPLEX[1,-c(1:2)])

gene_sets <- list(
  "GO:0070268 cornification (GOBP)" = GOBP_CORNIFICATION,

  "GO:0031424 keratinization (GOBP)" = GOBP_KERATINIZATION,

  "R-HSA-6809371 cornified envelope (REACTOME)" = REACTOME_CORNIFIED_ENVELOPLE,

  "Epidermal Differentiation Complex (1q21, Mischke et al.)" = EPIDERMAL_DIFFERENTIATION_COMPLEX
)
aliases <- list(LORICRIN = "LOR", LOR = "LORICRIN")

## ---- 2. Load DESeq2 results and define universe + DEGs ---------------------
res <- read.csv(deseq_csv, row.names = 1, stringsAsFactors = FALSE)

# Universe = genes actually tested. DESeq2 sets padj = NA for genes removed by
# independent filtering; those were never testable, so including them in the
# background inflates every p-value.
res_tested <- res[!is.na(res$padj), ]
universe   <- rownames(res_tested)
M          <- length(universe)

deg <- if (direction == "up") {
  rownames(res_tested)[res_tested$log2FoldChange >=  lfc_cut & res_tested$padj < fdr_cut]
} else {
  rownames(res_tested)[abs(res_tested$log2FoldChange) >= lfc_cut & res_tested$padj < fdr_cut]
}

cat(sprintf("Universe (padj != NA): %d/nDEGs (%s, |log2FC|>=%.2f, FDR<%.2f): %d/n  %s/n/n",
            M, direction, lfc_cut, fdr_cut, length(deg), paste(sort(deg), collapse = ", ")))


## ---- 3. Hypergeometric test ------------------------------------------------
# Map a set member onto the universe, trying aliases.
resolve <- function(g, universe) {
  if (g %in% universe) return(g)
  for (a in aliases[[g]]) if (!is.null(a) && a %in% universe) return(a)
  NA_character_
}

ora <- function(deg, gene_sets, universe) {
  M <- length(universe); k <- length(deg)
  out <- lapply(names(gene_sets), function(nm) {
    ann  <- unique(na.omit(vapply(gene_sets[[nm]], resolve, character(1), universe)))
    m    <- length(ann)                     # set members present in universe
    hits <- sort(intersect(deg, ann))
    q    <- length(hits)
    # P(X >= q); lower.tail=FALSE gives P(X > q-1). m = "white balls",
    # M - m = "black balls", k = balls drawn.
    p <- if (q > 0) phyper(q - 1, m, M - m, k, lower.tail = FALSE) else 1
    exp <- m * k / M
    data.frame(
      gene_set          = nm,
      genes_in_set      = length(gene_sets[[nm]]),
      genes_in_universe = m,
      n_DEG             = k,
      universe_size     = M,
      overlap           = q,
      expected_overlap  = round(exp, 4),
      fold_enrichment   = if (exp > 0) round(q / exp, 2) else NA_real_,
      p_hypergeometric  = p,
      overlap_genes     = paste(hits, collapse = "; "),
      stringsAsFactors  = FALSE)
  })
  out <- do.call(rbind, out)
  out <- out[order(out$p_hypergeometric), ]
  out$FDR_BH <- p.adjust(out$p_hypergeometric, method = "BH")
  rownames(out) <- NULL
  out
}

results <- ora(deg, gene_sets, universe)

print(results[, c("gene_set", "genes_in_universe", "overlap", "expected_overlap",
                  "fold_enrichment", "p_hypergeometric", "FDR_BH", "overlap_genes")],
      right = FALSE)

write.csv(results, paste0(out_prefix, "_results.csv"), row.names = FALSE)


## ---- 4. Sensitivity across thresholds --------------------------------------
# With only a handful of DEGs the p-values hinge on 2-3 genes, so check that
# the signal does not depend on where the fold-change cut lands.
thresholds <- list(
  "log2FC>=1 & FDR<0.1 (up)" = rownames(res_tested)[res_tested$log2FoldChange >= 1 & res_tested$padj < 0.1],
  "|log2FC|>=1 & FDR<0.1"    = rownames(res_tested)[abs(res_tested$log2FoldChange) >= 1 & res_tested$padj < 0.1],
  "FDR<0.1 & log2FC>0"       = rownames(res_tested)[res_tested$padj < 0.1 & res_tested$log2FoldChange > 0],
  "FDR<0.1 (any log2FC)"     = rownames(res_tested)[res_tested$padj < 0.1]
)
sens <- do.call(rbind, lapply(names(thresholds), function(lab) {
  d <- ora(thresholds[[lab]], gene_sets, universe); cbind(DEG_definition = lab, d)
}))
write.csv(sens, paste0(out_prefix, "_sensitivity.csv"), row.names = FALSE)

cat("/n-- Sensitivity across DEG definitions --/n")
print(sens[, c("DEG_definition", "gene_set", "overlap", "fold_enrichment", "FDR_BH")],
      right = FALSE)


## ---- 5. Per-gene detail for set members ------------------------------------
members <- unique(na.omit(vapply(unlist(gene_sets, use.names = FALSE),
                                resolve, character(1), universe)))
detail  <- res_tested[members, c("baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj")]
detail  <- detail[order(detail$padj), ]
detail$gene <- rownames(detail)
for (nm in names(gene_sets)) {
  ann <- na.omit(vapply(gene_sets[[nm]], resolve, character(1), universe))
  detail[[nm]] <- detail$gene %in% ann
}
write.csv(detail, paste0(out_prefix, "_gene_detail.csv"), row.names = FALSE)


## ---- 6. Bar plot -----------------------------------------------------------
bar_fill<-"#74BE5E"
v <- -log10(results$p_hypergeometric)
o <- order(v)
results <- results[o, ]; v <- v[o]
xmax <- ceiling(max(v))
pdf("ORA_barplot.pdf", width = 8.5, height = 3)
par(mar = c(4.2, 18, 2.8, 1.5),
    mgp = c(2.5, 0.6, 0), tcl = -0.3, xaxs = "i")
barplot(v, horiz = TRUE, names.arg = results$gene_set, las = 1, col = bar_fill,
        border = NA, xlim = c(0, xmax), axes = FALSE, cex.names = 0.85,
        space = 0.45,
        main = "Up DEGs Over-representation analysis")
axis(1, at = 0:xmax, cex.axis = 0.9)
mtext(expression(-log[10](italic(P) ~ value)), side = 1, line = 2.4)
# p = 0.05 reference. Delete this line for the clean look of the example figure.
abline(v = -log10(0.05), lty = 2, lwd = 3, col = "firebrick")
box()
dev.off()
