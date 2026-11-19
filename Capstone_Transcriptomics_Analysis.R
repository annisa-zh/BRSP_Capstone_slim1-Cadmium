# 00_PREPARASI DATA

install.packages("BiocManager")
BiocManager::install(c("GEOquery", "limma", "clusterProfiler", "org.At.tair.db"))
install.packages(c("ggplot2", "pheatmap"))

library(GEOquery)
library(limma)

gset <- getGEO("GSE138943", GSEMatrix = TRUE)[[1]]
ex <- exprs(gset)

dim(ex)
colnames(ex)

samples <- c("GSM4125057", "GSM4125058", "GSM4125059",
             "GSM4125054", "GSM4125055", "GSM4125056")
ex <- ex[, samples]

group <- factor(c("Control", "Control", "Control",
                  "Treatment", "Treatment", "Treatment"))

# 01_ANALISIS LIMMA
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
fit <- lmFit(ex, design)
kontras <- makeContrasts(Treatment - Control, levels = design)
fit2 <- eBayes(contrasts.fit(fit, kontras))
res <- topTable(fit2, number = Inf, adjust.method = "BH", sort.by = "P")
head(res)

# 02_TENTUKAN DEG

res$status <- "NO"

res$status[res$adj.P.Val < 0.05 & res$logFC >= 1] <- "UP"
res$status[res$adj.P.Val < 0.05 & res$logFC <= -1] <- "DOWN"

table(res$status)

deg <- res[res$status != "NO", ]
dim(deg)

dir.create("Hasil", showWarnings = FALSE)
write.csv(res, "Hasil/Tabel_hasil_semua_probe.csv")
write.csv(deg, "Hasil/Tabel_DEG.csv")

# 03_VOLCANO PLOT

library(ggplot2)

folder_hasil <- "C:/Users/Asus/Downloads/Hasil"
dir.create(folder_hasil, showWarnings = FALSE, recursive = TRUE)

volcano <- ggplot(res, aes(x = logFC, y = -log10(adj.P.Val), color = status)) +
  geom_point(alpha = 0.6, size = 1.2) +
  scale_color_manual(values = c(DOWN = "blue", NO = "grey70", UP = "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(title = "Volcano Plot: slim1-1 vs Col-0 pada 100 uM Cd (GSE138943)",
       x = "logFC",
       y = "-log10(adj. P value)") +
  theme_minimal()

volcano
ggsave(file.path(folder_hasil, "1_volcano_plot.png"),
       volcano, width = 7, height = 5, dpi = 300)


# 04_HEATMAP TOP 50

library(pheatmap)

top50 <- head(deg[order(deg$P.Value), ], 50)

mat <- ex[rownames(top50), ]

anotasi_kolom <- data.frame(Group = group, row.names = colnames(mat))

pheatmap(mat,
         scale = "row",
         clustering_distance_rows = "euclidean",
         clustering_method = "complete",
         annotation_col = anotasi_kolom,
         fontsize_row = 6,
         main = "Top 50 DEG: slim1-1 vs Col-0 (100 uM Cd)")

png(file.path(folder_hasil, "2_heatmap_top50.png"),
    width = 2000, height = 2400, res = 250)
pheatmap(mat,
         scale = "row",
         clustering_distance_rows = "euclidean",
         clustering_method = "complete",
         annotation_col = anotasi_kolom,
         fontsize_row = 6,
         main = "Top 50 DEG: slim1-1 vs Col-0 (100 uM Cd)")
dev.off()

# 05_ANOTASI GEN

anot <- fData(gset)
colnames(anot)

pola <- "AT[1-5CM]G[0-9]{5}"
skor <- sapply(anot, function(k) mean(grepl(pola, as.character(k), ignore.case = TRUE)))
skor
kolom_gen <- names(which.max(skor))
kolom_gen

ambil_locus <- function(x) {
  x <- as.character(x)
  m <- regmatches(x, regexpr(pola, x, ignore.case = TRUE))
  hasil <- rep(NA_character_, length(x))
  hasil[grepl(pola, x, ignore.case = TRUE)] <- toupper(m)
  hasil
}

res$Gene <- ambil_locus(anot[rownames(res), kolom_gen])

deg <- res[res$status != "NO", ]

head(res[, c("logFC", "adj.P.Val", "status", "Gene")])
sum(!is.na(res$Gene))                  # jumlah probe yang punya ID gen
length(unique(na.omit(deg$Gene)))      # jumlah gen unik di DEG
table(deg$status)                      # jumlah probe UP dan DOWN

write.csv(res, file.path(folder_hasil, "Tabel_hasil_semua_probe.csv"))
write.csv(deg, file.path(folder_hasil, "Tabel_DEG.csv"))

# 06_GO ENRICHMENT

library(clusterProfiler)
library(org.At.tair.db)
library(ggplot2)

folder_hasil <- "C:/Users/Asus/Downloads/Hasil"
dir.create(folder_hasil, showWarnings = FALSE, recursive = TRUE)

# Daftar gen unik dari DEG
gen_deg <- unique(na.omit(deg$Gene))
length(gen_deg)

# GO untuk tiga kategori
go_bp <- enrichGO(gene = gen_deg, OrgDb = org.At.tair.db, keyType = "TAIR",
                  ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05)
go_mf <- enrichGO(gene = gen_deg, OrgDb = org.At.tair.db, keyType = "TAIR",
                  ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.05)
go_cc <- enrichGO(gene = gen_deg, OrgDb = org.At.tair.db, keyType = "TAIR",
                  ont = "CC", pAdjustMethod = "BH", pvalueCutoff = 0.05)

# Jumlah term signifikan
nrow(as.data.frame(go_bp))
nrow(as.data.frame(go_mf))
nrow(as.data.frame(go_cc))

# Buat dotplot BP (15 term teratas)
dotplot_bp <- dotplot(go_bp, showCategory = 15) + ggtitle("GO Biological Process")
print(dotplot_bp)

# Simpan gambar
file_gambar <- file.path(folder_hasil, "3a_GO_BP_dotplot.png")
ggsave(filename = file_gambar, plot = dotplot_bp,
       width = 8, height = 7, dpi = 300)

# Simpan tabel
write.csv(as.data.frame(go_bp), file.path(folder_hasil, "Tabel_GO_BP.csv"))
write.csv(as.data.frame(go_mf), file.path(folder_hasil, "Tabel_GO_MF.csv"))
write.csv(as.data.frame(go_cc), file.path(folder_hasil, "Tabel_GO_CC.csv"))

# 07_KEGG PATHWAY

library(clusterProfiler)
library(ggplot2)

folder_hasil <- "C:/Users/Asus/Downloads/Hasil"
gen_deg <- unique(na.omit(deg$Gene))

kegg <- enrichKEGG(gene = gen_deg, organism = "ath", keyType = "kegg",
                   pAdjustMethod = "BH", pvalueCutoff = 0.05)

if (is.null(kegg) || nrow(as.data.frame(kegg)) == 0) {
  library(org.At.tair.db)
  konversi <- bitr(gen_deg, fromType = "TAIR", toType = "ENTREZID",
                   OrgDb = org.At.tair.db)
  kegg <- enrichKEGG(gene = unique(konversi$ENTREZID), organism = "ath",
                     keyType = "ncbi-geneid", pAdjustMethod = "BH",
                     pvalueCutoff = 0.05)
}

# Berapa jalur yang signifikan?
nrow(as.data.frame(kegg))


kegg@result$Description <- sub(" - Arabidopsis.*", "", kegg@result$Description)

# Lihat hasilnya
as.data.frame(kegg)[, c("Description", "GeneRatio", "Count", "p.adjust")]

# Dotplot dan simpan
dotplot_kegg <- dotplot(kegg, showCategory = 15) + ggtitle("KEGG Pathway")
print(dotplot_kegg)

file_kegg <- file.path(folder_hasil, "3b_KEGG_dotplot.png")
ggsave(filename = file_kegg, plot = dotplot_kegg, width = 8, height = 6, dpi = 300)

# Simpan tabel
write.csv(as.data.frame(kegg), file.path(folder_hasil, "Tabel_KEGG.csv"))

# Cek
file.exists(file_kegg)
list.files(folder_hasil)





#### RINGKASAN HASIL

up   <- deg[deg$status == "UP", ];   up   <- up[order(-up$logFC), ]
down <- deg[deg$status == "DOWN", ]; down <- down[order(down$logFC), ]

ringkasan <- data.frame(
  Keterangan = c("Jumlah probe dianalisis",
                 "Probe dengan adj.P.Val < 0,05",
                 "Probe DEG (adj.P < 0,05 dan |logFC| >= 1)",
                 "Probe UP",
                 "Probe DOWN",
                 "Probe yang punya ID gen",
                 "Gen unik DEG",
                 "Gen unik UP",
                 "Gen unik DOWN",
                 "Term GO BP signifikan",
                 "Term GO MF signifikan",
                 "Term GO CC signifikan",
                 "Jalur KEGG signifikan"),
  Nilai = c(nrow(res),
            sum(res$adj.P.Val < 0.05),
            nrow(deg),
            nrow(up),
            nrow(down),
            sum(!is.na(res$Gene)),
            length(unique(na.omit(deg$Gene))),
            length(unique(na.omit(up$Gene))),
            length(unique(na.omit(down$Gene))),
            nrow(as.data.frame(go_bp)),
            nrow(as.data.frame(go_mf)),
            nrow(as.data.frame(go_cc)),
            nrow(as.data.frame(kegg)))
)
print(ringkasan)

# Sepuluh probe UP dan DOWN dengan perubahan terbesar
head(up[,   c("Gene", "logFC", "P.Value", "adj.P.Val")], 10)
head(down[, c("Gene", "logFC", "P.Value", "adj.P.Val")], 10)

# Term GO dan jalur KEGG teratas
head(as.data.frame(go_bp)[, c("Description", "Count", "p.adjust")], 10)
as.data.frame(go_mf)[, c("Description", "Count", "p.adjust")]
as.data.frame(kegg)[, c("Description", "GeneRatio", "Count", "p.adjust")]

# Simpan ringkasan
write.csv(ringkasan, file.path(folder_hasil, "Ringkasan_angka.csv"), row.names = FALSE)
write.csv(head(up, 10),   file.path(folder_hasil, "Top10_UP.csv"))
write.csv(head(down, 10), file.path(folder_hasil, "Top10_DOWN.csv"))
