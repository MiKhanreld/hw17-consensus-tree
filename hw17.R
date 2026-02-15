library(ape)
library(purrr)
library(stringr)

setwd("C:/Users/Misha/Downloads/hw17")
file_path <- "table_with_frequencies.txt"

# 1) названия текстов (первая строка)
header_line <- readLines(file_path, n = 1, encoding = "UTF-8")
texts <- scan(text = header_line, what = character(), quiet = TRUE)

# 2) таблица
raw <- read.table(
  file_path,
  header = FALSE,
  skip = 1,
  quote = "\"",
  sep = "",
  stringsAsFactors = FALSE,
  check.names = FALSE
)
colnames(raw) <- c("word", texts)

# 3) матрица: тексты и слова
words <- raw$word
raw$word <- NULL

mx <- as.matrix(raw)
rownames(mx) <- make.unique(as.character(words))
X <- t(mx)  # тексты и слова

# 4) берем топ-500 самых частотных
top500 <- names(sort(colSums(X), decreasing = TRUE))[1:min(500, ncol(X))]
X500 <- X[, top500, drop = FALSE]

# 5) одно дерево: случайные 200 слов, z-score + Manhattan (Classic Delta)
get_tree <- function(Xm, mfw = 200) {
  cols <- sample(ncol(Xm), size = min(mfw, ncol(Xm)), replace = FALSE)
  X_sub <- Xm[, cols, drop = FALSE]
  Xz <- scale(X_sub)
  D <- dist(Xz, method = "manhattan")
  as.phylo(hclust(D, method = "average"))
}

set.seed(123)
trees <- map(1:100, ~get_tree(X500, mfw = 200))
cons <- consensus(trees, p = 0.5, rooted = FALSE)

# 6) цвета по авторам (часть до "_")
authors <- str_remove(cons$tip.label, "_.+$")
uniq_auth <- unique(authors)

pal <- c(
  "#1b9e77", "#d95f02", "#7570b3", "#e7298a",
  "#66a61e", "#e6ab02", "#a6761d", "#666666",
  "#1f78b4", "#b2df8a", "#fb9a99", "#cab2d6"
)
col_map <- setNames(pal[seq_along(uniq_auth)], uniq_auth)
tip_cols <- col_map[authors]

# 7) png: прямоугольное дерево, чтобы все влезало
png("consensus_tree.png", width = 1600, height = 1000, res = 150)

par(mar = c(3, 6, 4, 18)) # большой правый отступ под длинные подписи
plot.phylo(
  cons,
  type = "phylogram", # обычное дерево
  direction = "rightwards",
  use.edge.length = TRUE,
  font = 2,
  cex = 0.9,
  tip.color = tip_cols,
  edge.width = 1.2,
  label.offset = 0.01,
  no.margin = FALSE
)
title("Consensus tree (Classic Delta, 200 MFW, p=0.5)")

legend(
  "topleft",
  legend = names(col_map),
  col = col_map,
  pch = 16,
  bty = "n",
  cex = 0.9
)

dev.off()

cat(getwd(), "\n")
