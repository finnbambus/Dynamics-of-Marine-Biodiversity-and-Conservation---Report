# Pierre Noire macrobenthic time series (1977-2007): analysis behind the manuscript.
# Reproduces every figure and statistic reported in manuscript.typ.
# Input : "raw data/" (community, taxonomy and trait tables)
# Output: "figures/" (PNG files included by the manuscript)

suppressMessages({
  library(here)
  library(vegan)
  library(FD)
  library(strucchange)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(tidyr)
  library(viridisLite)
})
here::i_am("analysis.R")
dir.create(here("figures"), showWarnings = FALSE)

theme_paper <- theme_bw(base_size = 9) +
  theme(panel.border = element_blank(), axis.line = element_blank())
save_fig <- function(plot, name, width, height)
  ggsave(here("figures", name), plot, width = width, height = height, dpi = 300)
mako5 <- function(n = 5) viridisLite::mako(n, end = 0.9)

# ---- Data -------------------------------------------------------------------
abund_raw <- read.csv(here("raw data", "abundances-PN-77-07.csv"), sep = ";")
abund_raw[is.na(abund_raw)] <- 0
taxonomy  <- read.csv(here("raw data", "taxonomy-PN-77-07.csv"), sep = ";")
BDD       <- merge(taxonomy, abund_raw, by = "Species")   # all 381 species matched

year_cols <- grep("^X", names(abund_raw), value = TRUE)
dates     <- as.Date(sub("^X", "", year_cols), format = "%d.%m.%Y")
years     <- as.numeric(format(dates, "%Y"))
comm      <- t(as.matrix(abund_raw[, year_cols]))               # years x species
storage.mode(comm) <- "numeric"
colnames(comm) <- abund_raw$Species
S_full    <- as.matrix(BDD[, year_cols])                         # species x years (for levels and groups)
storage.mode(S_full) <- "numeric"
# SR and Shannon of a species x years matrix
idx_of <- function(M) data.frame(SR = colSums(M > 0), Shannon = apply(M, 2, function(x) vegan::diversity(x[x > 0], "shannon")))

# ---- Alpha diversity -------------------------------------------------------
SR      <- rowSums(comm > 0)
shannon <- vegan::diversity(comm, index = "shannon")
simpson <- vegan::diversity(comm, index = "simpson")
pielou  <- shannon / log(SR)
alpha   <- data.frame(SR, Shannon = shannon, Simpson = simpson, Pielou = pielou)

# ---- Dominance: k-dominance curves and top-5 species -----------------------
kdom <- do.call(rbind, lapply(seq_along(years), function(i) {
  x <- sort(comm[i, comm[i, ] > 0], decreasing = TRUE)
  data.frame(year = years[i], rank = seq_along(x), accumfreq = cumsum(x) / sum(x) * 100)
}))
kdom_plot <- ggplot(kdom, aes(rank, accumfreq, group = year, colour = year)) + geom_line() +
  xlab("Species rank") + ylab("Cumulative relative abundance (%)") +
  scale_x_continuous(trans = "log1p", breaks = c(1, 5, 10, 25, 50, 100)) +
  scale_colour_viridis_c(option = "mako", end = 0.9, name = "Year") +
  labs(tag = "a") + theme_paper +
  theme(plot.tag = element_text(face = "bold"), aspect.ratio = 1,
        legend.position = c(0.97, 0.05), legend.justification = c(1, 0), legend.direction = "horizontal",
        legend.background = element_rect(fill = alpha("white", 0.8), colour = NA),
        legend.key.width = unit(0.9, "cm"), legend.key.height = unit(0.3, "cm"), legend.title = element_text(size = 7), legend.text = element_text(size = 6.5))

top_species <- c("Ampelisca sarsi", "Spio decoratus", "Paradoneis armata",
                 "Ampelisca armoricana", "Ampelisca tenuicornis")
top5 <- do.call(rbind, lapply(seq_along(years), function(i) {
  x <- sort(comm[i, ], decreasing = TRUE)[1:5]
  data.frame(date = dates[i], species = names(x), abundance = as.numeric(x))
}))
total_abund <- data.frame(date = dates, total = rowSums(comm))
top5_sum    <- aggregate(abundance ~ date, top5, sum)
cmp         <- merge(total_abund, top5_sum, by = "date")
cor_top5    <- cor.test(cmp$total, cmp$abundance, method = "spearman", exact = FALSE)
cat(sprintf("Top-5 vs total: rho = %.3f, rho2 = %.2f, p = %.2e; mean share = %.1f%% (range %.1f-%.1f)\n",
            cor_top5$estimate, cor_top5$estimate^2, cor_top5$p.value, mean(100 * cmp$abundance / cmp$total),
            min(100 * cmp$abundance / cmp$total), max(100 * cmp$abundance / cmp$total)))
# Per year: leading species, share of the top 1, 2 and 5 species, and species needed to reach 95% of individuals
top_share <- t(apply(comm, 1, function(x) { x <- sort(x, decreasing = TRUE); 100 * cumsum(x)[c(1, 2, 5)] / sum(x) }))
n95 <- apply(comm, 1, function(x) which(cumsum(sort(x, decreasing = TRUE)) >= 0.95 * sum(x))[1])
print(data.frame(year = years, leading = colnames(comm)[apply(comm, 1, which.max)], top1 = round(top_share[, 1]),
                 top2 = round(top_share[, 2]), top5 = round(top_share[, 3], 1), n95 = n95), row.names = FALSE)
cat("Years needing more than 20 species to reach 95%:", years[n95 > 20], "\n")

# species order and colours follow the 1977 bar (largest first); "Other" stays at the bottom
top5$grp <- ifelse(top5$species %in% top_species, top5$species, "Other")
first_bar <- top5[top5$date == min(top5$date) & top5$grp != "Other", ]
first_order <- as.character(first_bar$grp[order(-first_bar$abundance)])
lvl <- c(first_order, setdiff(top_species, first_order), "Other")
top5$grp <- factor(top5$grp, levels = lvl)
bar_cols <- c(setNames(rev(mako5()), lvl[lvl != "Other"]), Other = "grey70")
stack <- top5 %>% group_by(date) %>% arrange(desc(grp == "Other"), abundance, .by_group = TRUE) %>%
  mutate(ymax = cumsum(abundance), ymin = ymax - abundance) %>% ungroup()

top5_plot <- ggplot() +
  geom_line(data = total_abund, aes(date, total), colour = "black", linewidth = 0.5) +
  geom_point(data = total_abund, aes(date, total), colour = "black", size = 0.9) +
  geom_rect(data = stack, aes(xmin = date - 165, xmax = date + 165, ymin = ymin, ymax = ymax, fill = grp)) +
  scale_fill_manual(values = bar_cols, name = NULL) +
  scale_y_continuous(labels = function(x) x / 1e4, expand = expansion(mult = c(0.02, 0.3))) +
  xlab("Sampling year") + ylab(expression(paste("Abundance (", 10^4, " ind.)"))) +
  labs(tag = "b") + theme_paper +
  theme(plot.tag = element_text(face = "bold"), legend.position = c(0.02, 0.98), legend.justification = c(0, 1), legend.direction = "horizontal",
        legend.background = element_rect(fill = alpha("white", 0.8), colour = NA),
        legend.text = element_text(size = 6.5), legend.key.size = unit(0.28, "cm"),
        legend.spacing.x = unit(0.1, "cm")) +
  guides(fill = guide_legend(ncol = 2, byrow = TRUE))

# k-dominance (square) and abundance (rest of the width) side by side
png(here("figures", "dominance_structure.png"), width = 7.2, height = 3.0, units = "in", res = 300)
grid::grid.newpage()
grid::pushViewport(grid::viewport(layout = grid::grid.layout(1, 2, widths = grid::unit(c(3.0, 1), c("in", "null")))))
print(kdom_plot, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
print(top5_plot, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
invisible(dev.off())

# ---- Taxonomic distinctness (Delta+) ---------------------------------------
tax_dist <- taxa2dist(taxonomy, varstep = TRUE)
dplus <- as.numeric(taxondive(comm, tax_dist)$Dplus)

# Per-year values and summaries of the alpha indices quoted in the text
cat(sprintf("Data: %d annual samples (%d-%d), %d species\n", nrow(comm), min(years), max(years), ncol(comm)))
print(data.frame(year = years, SR, Shannon = round(shannon, 2), Dplus = round(dplus, 1), Abundance = rowSums(comm)), row.names = FALSE)
cat(sprintf("SR: range %d-%d, mean %.0f, SD %.0f\n", min(SR), max(SR), mean(SR), sd(SR)))
cat(sprintf("Delta+: range %.1f-%.1f, mean %.1f; linear trend over time R2 = %.3f\n",
            min(dplus), max(dplus), mean(dplus), summary(lm(dplus ~ years))$r.squared))

# ---- Alpha-diversity panel: 3 x 2 grid (5 time series + correlation table) -
# hl = TRUE marks the indices carried into the main discussion: tinted frame and accent-coloured series
hl_col <- mako5()[3]
hl_fill <- rgb(t(col2rgb(hl_col) * 0.08 + 255 * 0.92), maxColorValue = 255)   # opaque card colour (8% tint)
ts_panel <- function(y, ylab, tag, hl = FALSE)
  ggplot(data.frame(date = dates, y = y), aes(date, y)) +
    geom_line(linewidth = 0.3, colour = if (hl) hl_col else "black") +
    geom_point(size = 1.2, colour = if (hl) hl_col else "black") +
    labs(x = "Sampling year", y = ylab, tag = tag) + theme_paper +
    theme(plot.tag = element_text(face = "bold"), plot.margin = margin(8, 10, 6, 8),
          plot.background = element_blank(),            # frame for hl panels is drawn in save_grid (rounded)
          panel.background = element_rect(fill = if (hl) hl_fill else "white", colour = NA))
# Spearman correlation heatmap: numbers are rho^2 (shared variance), colour is -log10(p) in mako (darker = smaller p),
# grey = not significant (p >= 0.05), white diagonal
# Width (pt) of the widest row label; the left margin is padded so every matrix gets the same tile size
label_w <- function(x) {
  pdf(NULL); on.exit(dev.off())                                   # off-screen device for measuring
  max(sapply(x, function(s)
    grid::convertWidth(grid::grobWidth(grid::textGrob(s, gp = grid::gpar(fontsize = 7))), "pt", valueOnly = TRUE)))
}
cor_heatmap <- function(df, labels, tag) {
  nm   <- names(df)
  pad  <- max(0, label_w("Shannon") - label_w(labels))
  cor_m <- cor(df, method = "spearman")
  p_m  <- outer(nm, nm, Vectorize(function(r, c)
    if (r == c) NA else cor.test(df[[r]], df[[c]], method = "spearman", exact = FALSE)$p.value))
  dimnames(p_m) <- dimnames(cor_m)
  print(round(cor_m, 2)); print(round(cor_m^2, 2)); print(signif(p_m, 2))
  d <- expand.grid(row = nm, col = nm, stringsAsFactors = FALSE)
  d$rho  <- mapply(function(r, c) cor_m[r, c], d$row, d$col)
  d$p    <- mapply(function(r, c) p_m[r, c], d$row, d$col)
  d$diag <- d$row == d$col
  d$logp <- ifelse(!is.na(d$p) & d$p < 0.05, -log10(d$p), NA)
  d$row  <- factor(d$row, levels = rev(nm))
  d$col  <- factor(d$col, levels = nm)
  ggplot(d, aes(col, row, fill = logp)) +
    geom_tile(colour = "white", linewidth = 0.8) +
    geom_tile(data = subset(d, diag), fill = "white", colour = "white", linewidth = 0.8) +
    geom_text(aes(label = sprintf("%.2f", rho^2), colour = !is.na(logp) & logp > 12), size = 3) +
    scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black"), guide = "none") +
    scale_fill_viridis_c(option = "mako", begin = 0.05, end = 0.95, direction = -1, na.value = "grey85",
                         limits = c(-log10(0.05), 20), oob = scales::squish, guide = "none") +
    scale_x_discrete(position = "top", labels = labels, expand = c(0, 0)) +
    scale_y_discrete(labels = rev(labels), expand = c(0, 0)) +
    labs(x = NULL, y = NULL, tag = tag) + theme_paper +
    theme(plot.tag = element_text(face = "bold"), axis.ticks = element_blank(),
          axis.text = element_text(size = 7),
          plot.margin = margin(5.5, 28, 22, 5.5 + pad))
}
# One grid cell: panel letter outside, plot (optionally in a rounded highlight card) to its right/below
draw_cell <- function(p, tag, row, col, hl = FALSE, plot_h = 0.90, card_h = 0.98, w = 0.925, card_w = w) {
  grid::pushViewport(grid::viewport(layout.pos.row = row, layout.pos.col = col))
  grid::grid.text(tag, x = grid::unit(0.012, "npc"), y = grid::unit(0.985, "npc"), just = c("left", "top"),
                  gp = grid::gpar(fontface = "bold", fontsize = 11))
  if (hl) {                                                                  # card reaches higher than the plot (room at the top)
    grid::pushViewport(grid::viewport(x = 0.06, y = 0.01, width = card_w, height = card_h, just = c("left", "bottom")))
    grid::grid.roundrect(r = grid::unit(0.07, "snpc"), gp = grid::gpar(fill = hl_fill, col = hl_col, lwd = 1.6))
    grid::popViewport()
  }
  grid::pushViewport(grid::viewport(x = 0.06, y = 0.01, width = w, height = plot_h, just = c("left", "bottom")))
  print(p + labs(tag = NULL), newpage = FALSE)
  grid::popViewport(2)
}
save_grid <- function(panels, name, height, hl = integer()) {
  png(here("figures", name), width = 7.2, height = height, units = "in", res = 300)
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(3, 2)))
  for (i in seq_along(panels))
    draw_cell(panels[[i]], letters[i], (i - 1) %/% 2 + 1, (i - 1) %% 2 + 1, hl = i %in% hl)
  invisible(dev.off())
}
alpha_panels <- list(
  ts_panel(SR, "Species richness", "a", hl = TRUE),
  ts_panel(shannon, "Shannon H'", "b", hl = TRUE),
  ts_panel(simpson, "Simpson 1 - D", "c"),
  ts_panel(pielou, "Pielou J", "d"),
  ts_panel(dplus, expression(paste("Tax. distinctness ", Delta^+{})), "e"),
  cor_heatmap(cbind(alpha, Dplus = dplus, Abundance = rowSums(comm)),
               c("SR", "Shannon", "Simpson", "Pielou", "\u0394+", "Abund."), "f"))
save_grid(alpha_panels, "alpha_diversity.png", 4.46, hl = 1:2)

# ---- Annex: Hill numbers of order q = 0, 1, 2 (effective numbers of species) --
p_rel <- comm / rowSums(comm)
hill  <- data.frame(date = dates, q0 = SR, q1 = exp(shannon), q2 = 1 / rowSums(p_rel^2))
hill_long <- pivot_longer(hill, -date, names_to = "order", values_to = "D")
hill_lab  <- c(q0 = "q = 0 (richness)", q1 = "q = 1 (exp Shannon)", q2 = "q = 2 (inverse Simpson)")
hill_long$order <- factor(hill_lab[hill_long$order], levels = hill_lab)
hill_plot <- ggplot(hill_long, aes(date, D, colour = order, shape = order)) +
  geom_line(linewidth = 0.4) + geom_point(size = 1.8) +
  scale_colour_manual(values = setNames(mako5(3), hill_lab)) +
  xlab("Sampling year") + ylab("Hill number (effective number of species)") +
  theme_paper + theme(legend.position = "top", legend.title = element_blank())
save_fig(hill_plot, "hill_numbers.png", 5.4, 3.6)

# ---- Beta diversity: Bray-Curtis and Jaccard ordinations ----------------------
X <- log1p(comm)
bray <- vegdist(X, method = "bray")
set.seed(1)
nmds_bray <- metaMDS(bray, k = 2, autotransform = FALSE, trace = FALSE)
cat("Bray-Curtis stress:", round(nmds_bray$stress, 3), "\n")
sc <- as.data.frame(scores(nmds_bray, display = "sites"))
sc$date <- dates

hc   <- hclust(bray, method = "average")

# Phases: UPGMA tree cut into five groups (checked on the nMDS and dendrogram)
clus <- cutree(hc, k = 5)
phase_levels <- c("Pre-spill", "Post-spill crash", "Post-spill remnant", "Recovery", "Millennium decline")
phase_cols   <- setNames(mako5(), phase_levels)
sc$phase <- factor(phase_levels[clus], levels = phase_levels)

# Check: breakpoints in the mean of SR and of Shannon (BIC, minimum segment 4 years; confint = 95% interval of the break year).
# A phase boundary is attributed to a series when the last year of the old phase lies inside that interval.
# The spill cannot be tested because the series starts one year before it.
bp_tab <- lapply(list(SR = SR, `H'` = shannon), function(y) {
  bp <- breakpoints(y ~ 1, h = 4)
  k  <- which.min(summary(bp)$RSS["BIC", ]) - 1
  if (k == 0) return(NULL)
  ci <- matrix(years[pmax(1, confint(bp, breaks = k)$confint)], nrow = k)  # lower, break, upper (lower capped at 1977)
  colnames(ci) <- c("lower", "break", "upper"); ci
})
for (s in names(bp_tab)) { cat(s, "breakpoints (last year of old regime):\n"); print(bp_tab[[s]]) }
prev_year <- sapply(seq_along(phase_levels), function(g) if (g == 1) NA else max(years[clus == g - 1]))
driver <- sapply(prev_year, function(y) {
  if (is.na(y)) return("")
  hit <- names(bp_tab)[sapply(bp_tab, function(ci) !is.null(ci) && any(y >= ci[, "lower"] & y <= ci[, "upper"]))]
  if (length(hit)) paste(hit, collapse = " + ") else "none" })
driver[2] <- "n/a"                                                 # spill: not testable
cat("Series with a breakpoint at each phase start:", paste(driver, collapse = " | "), "\n")
# Mean richness, Shannon diversity, abundance and leading taxa per phase
for (g in seq_along(phase_levels)) {
  m <- colMeans(comm[clus == g, , drop = FALSE]); t4 <- sort(m, decreasing = TRUE)[1:4]
  cat(sprintf("%s: mean SR %.0f, Shannon %.2f, abundance %.0f; leading taxa %s\n", phase_levels[g],
              mean(SR[clus == g]), mean(shannon[clus == g]), sum(m),
              paste0(names(t4), " (", round(100 * t4 / sum(m)), "%)", collapse = ", ")))
}

starts <- sc[c(TRUE, head(sc$phase, -1) != tail(sc$phase, -1)), ]
offsets <- data.frame(year = c("1977", "1978", "1982", "1993", "2005"),
                      dx = c(-0.02, 0, 0.035, -0.035, -0.035),
                      dy = c(-0.035, -0.035, -0.035, 0.035, 0))
starts$year <- format(starts$date, "%Y")
starts <- merge(starts, offsets, by = "year")
spill <- sc[format(sc$date, "%Y") == "1978", ]

bray_plot <- ggplot(sc, aes(NMDS1, NMDS2)) +
  geom_path(arrow = arrow(type = "closed", length = unit(0.03, "npc")), colour = "grey50") +
  geom_point(aes(colour = phase), size = 3) +
  geom_text(data = subset(starts, year != "1978"), aes(x = NMDS1 + dx, y = NMDS2 + dy, label = format(date, "%Y")), size = 3) +
  geom_segment(data = spill, aes(x = NMDS1, xend = NMDS1, y = NMDS2 - 0.012, yend = NMDS2 - 0.10),
               colour = phase_cols[["Post-spill crash"]], linewidth = 0.6) +
  geom_label(data = spill, aes(x = NMDS1, y = NMDS2 - 0.118, label = "Amoco Cadiz oil spill (1978)"),
             fill = phase_cols[["Post-spill crash"]], colour = "white", fontface = "bold", size = 3) +
  expand_limits(y = -0.44, x = 0.42) +
  scale_colour_manual(values = phase_cols) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5, size = 3,
           label = paste0("Stress = ", round(nmds_bray$stress, 3))) +
  xlab("nMDS1") + ylab("nMDS2") + theme_paper + theme(legend.position = "none", aspect.ratio = 1)
bray_plot <- bray_plot + labs(tag = "a") + theme(plot.tag = element_text(face = "bold", size = 11))

# Presence/absence ordinations (annex): Jaccard and Sorensen, square panels, half the width each
pa_plot <- function(d, name, tag) {
  set.seed(1)
  nm <- metaMDS(d, k = 2, autotransform = FALSE, trace = FALSE)
  cat(name, "stress:", round(nm$stress, 3), "\n")
  sc_pa <- as.data.frame(scores(nm, display = "sites")); sc_pa$date <- dates
  ggplot(sc_pa, aes(NMDS1, NMDS2)) +
    geom_path(arrow = arrow(type = "closed", length = unit(0.03, "npc")), colour = "grey50") +
    geom_point(size = 1.8) +
    geom_text_repel(aes(label = format(date, "%Y")), size = 2.5, max.overlaps = Inf) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5, size = 3,
             label = paste0("Stress = ", round(nm$stress, 3))) +
    labs(title = name, tag = tag, x = "nMDS1", y = "nMDS2") + theme_paper +
    theme(aspect.ratio = 1, plot.tag = element_text(face = "bold", size = 11),
          plot.title = element_text(face = "bold", size = 9, hjust = 0.5))
}
png(here("figures", "nmds_presence_absence.png"), width = 7.2, height = 3.8, units = "in", res = 300)
grid::grid.newpage()
grid::pushViewport(grid::viewport(layout = grid::grid.layout(1, 2)))
print(pa_plot(vegdist(comm, method = "jaccard", binary = TRUE), "Jaccard", "a"),
      vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
print(pa_plot(vegdist(comm, method = "bray", binary = TRUE), "S\u00f8rensen", "b"),
      vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
invisible(dev.off())

# UPGMA dendrogram coloured by phase (branches and leaves)
dend <- as.dendrogram(hc)
year_lab <- setNames(format(dates, "%Y"), rownames(comm))
color_branches <- function(node) {
  if (is.leaf(node)) {
    orig <- attr(node, "label"); g <- clus[[orig]]
    attr(node, "label") <- year_lab[[orig]]
    attr(node, "edgePar") <- list(col = phase_cols[[phase_levels[g]]])
    attr(node, "cluster") <- g
    return(node)
  }
  node[[1]] <- color_branches(node[[1]]); node[[2]] <- color_branches(node[[2]])
  g1 <- attr(node[[1]], "cluster"); g2 <- attr(node[[2]], "cluster")
  if (!is.na(g1) && !is.na(g2) && g1 == g2) {
    attr(node, "cluster") <- g1
    attr(node, "edgePar") <- list(col = phase_cols[[phase_levels[g1]]])
  } else {
    attr(node, "cluster") <- NA
    attr(node, "edgePar") <- list(col = "grey50")
  }
  node
}
dend <- color_branches(dend)
# Phase summary table: years and total-abundance range per phase
tot_ab <- rowSums(comm)
phase_tab <- do.call(rbind, lapply(seq_along(phase_levels), function(g) {
  i <- which(clus == g)
  yr <- range(years[i]); ab <- range(tot_ab[i])
  data.frame(Period = phase_levels[g],
             Years = if (yr[1] == yr[2]) as.character(yr[1]) else paste0(yr[1], "\u2013", yr[2]),
             Abundance = if (ab[1] == ab[2]) format(ab[1]) else paste0(format(ab[1]), "\u2013", format(ab[2])),
             Break = driver[g])
}))
print(phase_tab)
tab_long <- do.call(rbind, lapply(seq_len(nrow(phase_tab)), function(r)
  data.frame(row = r, col = c("Period", "Years", "Abundance", "Break in"),
             label = c(phase_tab$Period[r], phase_tab$Years[r], phase_tab$Abundance[r], phase_tab$Break[r]),
             colour = phase_cols[[r]])))
col_x <- c(Period = 0, Years = 0.34, Abundance = 0.60, `Break in` = 0.87)
tab_long$x <- col_x[tab_long$col]
tab_long$colour[tab_long$col != "Period"] <- "black"
tab_plot <- ggplot() +
  geom_text(data = tab_long, aes(x, -row, label = label, colour = I(colour)), hjust = 0, size = 2.9) +
  geom_text(data = data.frame(x = col_x, label = names(col_x)),
            aes(x, 0.2, label = label), hjust = 0, size = 2.9, fontface = "bold") +
  annotate("segment", x = -0.02, xend = 1.02, y = -0.35, yend = -0.35, linewidth = 0.3) +
  scale_x_continuous(limits = c(-0.03, 1.03), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-5.6, 0.7), expand = c(0, 0)) + theme_void()

# UPGMA dendrogram coloured by phase, drawn as a ggplot (segments from the coloured dendrogram)
dend_seg <- list(); dend_lab <- list(); leaf_i <- 0
walk <- function(node) {                          # returns c(x, height) of the node
  if (is.leaf(node)) {
    leaf_i <<- leaf_i + 1
    dend_lab[[length(dend_lab) + 1]] <<- data.frame(x = leaf_i, label = attr(node, "label"),
                                                    colour = attr(node, "edgePar")$col)
    return(c(leaf_i, 0))
  }
  h <- attr(node, "height")
  kids <- lapply(seq_len(2), function(k) walk(node[[k]]))
  cols <- sapply(1:2, function(k) { e <- attr(node[[k]], "edgePar"); if (is.null(e)) "grey50" else e$col })
  xm <- mean(c(kids[[1]][1], kids[[2]][1]))
  for (k in 1:2) {
    dend_seg[[length(dend_seg) + 1]] <<- data.frame(x = kids[[k]][1], xend = kids[[k]][1], y = h, yend = kids[[k]][2], colour = cols[k])
    dend_seg[[length(dend_seg) + 1]] <<- data.frame(x = xm, xend = kids[[k]][1], y = h, yend = h, colour = cols[k])
  }
  c(xm, h)
}
invisible(walk(dend))
dend_seg <- do.call(rbind, dend_seg); dend_lab <- do.call(rbind, dend_lab)
dend_plot <- ggplot() +
  geom_segment(data = dend_seg, aes(x, y, xend = xend, yend = yend, colour = I(colour)), linewidth = 0.5, lineend = "round") +
  geom_text(data = dend_lab, aes(x, -0.01, label = label, colour = I(colour)), angle = 90, hjust = 1, size = 2.2) +
  scale_y_continuous(limits = c(-0.095, 0.4), breaks = seq(0, 0.4, 0.1), expand = c(0, 0)) +
  scale_x_continuous(expand = expansion(add = 0.8)) +
  ylab("Bray-Curtis dissimilarity") + labs(tag = "c") + theme_paper +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.line.y = element_line(), panel.grid = element_blank(), plot.tag = element_text(face = "bold", size = 11))

# Beta-diversity figure: nMDS (left); phase table over UPGMA dendrogram (right)
png(here("figures", "beta_diversity.png"), width = 7.2, height = 3.6, units = "in", res = 300)
grid::grid.newpage()
grid::pushViewport(grid::viewport(layout = grid::grid.layout(2, 2, widths = grid::unit(c(1, 1), "null"),
                                                             heights = grid::unit(c(1.25, 1), "null"))))
print(bray_plot, vp = grid::viewport(layout.pos.row = 1:2, layout.pos.col = 1))
print(tab_plot + labs(tag = "b") + theme(plot.tag = element_text(face = "bold", size = 11), plot.margin = margin(5.5, 8, 0, 5.5)),
      vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
print(dend_plot, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 2))
invisible(dev.off())

# ---- Functional diversity ---------------------------------------------------
trait_tab <- read.csv(here("raw data", "FC-PN-77-07.csv"), sep = ";", row.names = 1)
ab_fd     <- read.csv(here("raw data", "Abund-PN-77-07.csv"), sep = ";", row.names = 1)
trait_blocks <- c(5, 3, 6, 2, 3, 6, 6, 5)                           # modalities per trait
cat(sprintf("Traits: %d taxa, %d traits, %d modalities (%d of %d species have no trait entry)\n",
            nrow(trait_tab), length(trait_blocks), sum(trait_blocks), ncol(comm) - nrow(trait_tab), ncol(comm)))
trait     <- prep.fuzzy.var(trait_tab, trait_blocks, row.w = rep(1, nrow(trait_tab)))
fd <- suppressWarnings(dbFD(trait, ab_fd, m = "min", w.abun = TRUE, stand.FRic = FALSE,
                            stand.x = FALSE, scale.RaoQ = FALSE,
                            calc.FDiv = TRUE, calc.CWM = FALSE))
# Spearman correlations of the main functional indices with the alpha indices
alpha_ref <- cbind(alpha, Dplus = dplus)
fd_alpha <- do.call(rbind, lapply(c("RaoQ", "FRic", "FEve"), function(f) do.call(rbind, lapply(names(alpha_ref), function(a) {
  tt <- cor.test(fd[[f]], alpha_ref[[a]], method = "spearman", exact = FALSE)
  data.frame(fd = f, alpha = a, rho = round(unname(tt$estimate), 2), rho2 = round(unname(tt$estimate)^2, 2), p = signif(tt$p.value, 2))
}))))
print(fd_alpha)
# Per-year values and ranges of the main functional indices
print(data.frame(year = years, RaoQ = round(fd$RaoQ, 2), FRic = round(fd$FRic, 2), FEve = round(fd$FEve, 3)), row.names = FALSE)
cat(sprintf("Rao's Q: %.2f (%d) to %.2f (%d); FRic: %.2f-%.2f, change 1977-1978 %.0f%%; FEve: %.2f-%.2f\n",
            min(fd$RaoQ), years[which.min(fd$RaoQ)], max(fd$RaoQ), years[which.max(fd$RaoQ)],
            min(fd$FRic), max(fd$FRic), 100 * (fd$FRic[2] / fd$FRic[1] - 1), min(fd$FEve), max(fd$FEve)))

# Functional-diversity panel: 3 x 2 grid (5 time series + correlation table)
fd_df <- data.frame(RaoQ = fd$RaoQ, FRic = fd$FRic, FEve = fd$FEve, FDiv = fd$FDiv, FDis = fd$FDis)
fd_cor_df <- cbind(fd_df, Abundance = rowSums(comm))
fd_panels <- list(
  ts_panel(fd_df$RaoQ, "Rao's Q", "a", hl = TRUE),
  ts_panel(fd_df$FRic, "Functional richness", "b", hl = TRUE),
  ts_panel(fd_df$FEve, "Functional evenness", "c", hl = TRUE),
  ts_panel(fd_df$FDiv, "Functional divergence", "d"),
  ts_panel(fd_df$FDis, "Functional dispersion", "e"),
  cor_heatmap(fd_cor_df, c("Rao Q", "FRic", "FEve", "FDiv", "FDis", "Abund."), "f"))
save_grid(fd_panels, "functional_diversity.png", 4.46, hl = 1:3)

# CWM of six traits, dominant species removed
blocks <- list(size = paste0("A.", 1:5), lifespan = paste0("B.", 1:3),
               feeding = paste0("C.", 1:6), development = paste0("D.", 1:3),
               movement = paste0("E.", 1:6), tolerance = paste0("G.", 1:5))
titles <- c(size = "Maximum body size", lifespan = "Life span", feeding = "Feeding mode",
            development = "Development type", movement = "Movement type",
            tolerance = "Tolerance to disturbance (AMBI)")
cat_labels <- list(
  size = c("Very small (<1 cm)", "Small (1-2 cm)", "Small-medium (3-10 cm)", "Medium (11-20 cm)", "Medium-large (>20 cm)"),
  lifespan = c("Short (<2 y)", "Medium (2-5 y)", "Long (>5 y)"),
  feeding = c("Filter-feeder", "Surface deposit", "Sub-surface deposit", "Predator/scavenger", "Omnivore", "Grazer"),
  development = c("Direct development", "Planktotrophic larva", "Lecithotrophic larva"),
  movement = c("Burrower", "Crawler", "Swimmer", "Walker", "Sessile", "Tube dweller"),
  tolerance = c("Sensitive", "Indifferent", "Tolerant", "2nd-order opportunist", "1st-order opportunist"))
trait_fuzzy <- function(cols) {
  m <- as.matrix(trait_tab[, cols]); m[is.na(m)] <- 0
  rs <- rowSums(m); rs[rs == 0] <- NA; sweep(m, 1, rs, "/")
}
top_ab <- c("Ampelisca_sarsi", "Spio_decorata", "Paradoneis_armata", "Ampelisca_armoricana", "Ampelisca_tenuicornis")
keep <- setdiff(colnames(ab_fd), top_ab)
P <- as.matrix(ab_fd[, keep]); P[is.na(P)] <- 0; P <- P / rowSums(P)
cwm_long <- bind_rows(lapply(names(blocks), function(t) {
  F <- trait_fuzzy(blocks[[t]])[keep, , drop = FALSE]; F[is.na(F)] <- 0
  m <- P %*% F; colnames(m) <- cat_labels[[t]]
  df <- as.data.frame(m); df$year <- years; df$trait <- titles[[t]]
  pivot_longer(df, cols = -c(year, trait), names_to = "category", values_to = "cwm") %>%
    mutate(category = factor(category, levels = cat_labels[[t]]),
           colour = setNames(mako5(length(cat_labels[[t]])), cat_labels[[t]])[as.character(category)])
}))
cwm_long$trait <- factor(cwm_long$trait, levels = unname(titles))
cwm_long$key <- paste(cwm_long$trait, cwm_long$category)
# Year and value of the maximum share of each trait category, and the full series for sensitive taxa
print(as.data.frame(cwm_long %>% group_by(trait, category) %>%
  summarise(max_year = year[which.max(cwm)], max = round(max(cwm), 2), .groups = "drop")))
with(subset(cwm_long, category == "Sensitive"), cat("Sensitive share:", paste0(year, " ", round(cwm, 2), collapse = ", "), "\n"))
cwm_ends <- cwm_long %>% group_by(trait, category) %>% filter(year == max(year)) %>% ungroup()
cwm_plot <- ggplot(cwm_long, aes(year, cwm, group = key, colour = I(colour))) +
  geom_line(linewidth = 0.5) + geom_point(size = 0.9) +
  geom_text_repel(data = cwm_ends, aes(label = category), hjust = 0, nudge_x = 0.8, size = 2.4,
                  direction = "y", segment.size = 0.2, box.padding = 0.08, min.segment.length = 0.3, seed = 1) +
  facet_wrap(~trait, ncol = 2, scales = "free_y") +
  scale_x_continuous(breaks = seq(1980, 2000, 10), expand = expansion(mult = c(0.02, 0.62))) +
  xlab("Sampling year") + ylab("Community-weighted mean (share of community)") + theme_paper +
  theme(legend.position = "none", strip.background = element_blank(), strip.text = element_text(face = "bold"))
save_fig(cwm_plot, "cwm_traits.png", 7.28, 5.6)

# ---- Taxonomic sufficiency: higher levels against species -------------------
levels_tax <- c("Genus", "Family", "Order", "Class", "Phylum")
species_idx <- idx_of(S_full)
level_matrix <- function(lev) {                                     # taxa x years, abundances summed per taxon
  agg <- aggregate(BDD[, year_cols], by = list(taxon = BDD[[lev]]), FUN = sum)
  m <- as.matrix(agg[, year_cols]); storage.mode(m) <- "numeric"; m
}
level_all <- lapply(setNames(levels_tax, levels_tax), level_matrix)
n_taxa <- sapply(level_all, nrow)
print(data.frame(level = levels_tax, taxa = n_taxa, pct_of_species = round(100 * n_taxa / ncol(comm))), row.names = FALSE)
suff <- bind_rows(lapply(levels_tax, function(l) {
  L <- idx_of(level_all[[l]])
  bind_rows(lapply(c("SR", "Shannon"), function(ix) {
    fit <- lm(L[[ix]] ~ species_idx[[ix]]); sm <- summary(fit)
    data.frame(level = l, index = ix,
               r = round(sqrt(sm$r.squared) * sign(coef(fit)[2]), 3),
               R2 = round(sm$r.squared, 3), p = signif(coef(sm)[2, 4], 3))
  }))
}))
print(suff)
suff_long <- bind_rows(lapply(levels_tax, function(l)
  data.frame(level = l, species_SR = species_idx$SR, SR = idx_of(level_all[[l]])$SR)))
# Bray-Curtis ordinations at each level (stress in panel titles) and Mantel tests
ord_all <- list(); bc_all <- list(); stress_all <- c()
for (lev in c("Species", levels_tax)) {
  X <- if (lev == "Species") comm else t(level_all[[lev]])
  d <- vegdist(log1p(X), method = "bray"); set.seed(1)
  nm <- metaMDS(d, k = 2, autotransform = FALSE, trace = FALSE)
  s <- as.data.frame(scores(nm, display = "sites")); s$year <- years; s$level <- lev
  ord_all[[lev]] <- s; bc_all[[lev]] <- d; stress_all[lev] <- nm$stress
}
mantel_tab <- do.call(rbind, lapply(levels_tax, function(l) {
  mt <- mantel(bc_all$Species, bc_all[[l]], method = "pearson", permutations = 999)
  data.frame(level = l, r = round(mt$statistic, 3), R2 = round(mt$statistic^2, 3), p = signif(mt$signif, 3))
}))
print(mantel_tab); print(round(stress_all, 3))

# ---- Sufficiency figure: two 3 x 2 grids side by side ---------------------------
# left: SR at each level vs species-level SR (square panels) + Shannon results table
# right: Bray-Curtis nMDS at each level (Class ... Species) + Mantel results table
level_cols <- setNames(mako5(), c("Phylum", "Class", "Order", "Family", "Genus"))
tag_i <- 0
next_tag <- function() { tag_i <<- tag_i + 1; letters[tag_i] }
tag_theme <- theme(plot.tag = element_text(face = "bold", size = 11),
                   plot.title = element_text(face = "bold", size = 8.5, hjust = 0.5))

focus_level <- "Family"                                              # level highlighted in the figure
suff_panel <- function(l) {
  d  <- subset(suff_long, level == l)
  ggplot(d, aes(species_SR, SR)) +
    geom_point(size = 1.4, colour = level_cols[[l]]) +
    geom_smooth(method = "lm", se = TRUE, level = 0.95, colour = "grey40", fill = "grey60", alpha = 0.25, linewidth = 0.5, formula = y ~ x) +
    labs(x = "SR (species level)", y = paste("SR at", tolower(l), "level"), title = l, tag = next_tag()) +
    theme_paper + tag_theme + theme(axis.title = element_text(size = 7.5),
                                    axis.text = element_text(size = 7),
                                    plot.background = element_blank(),
                                    panel.background = element_rect(fill = if (l == focus_level) hl_fill else "white", colour = NA))
}

res_table <- function(df, title, tag) {                           # df: Level, r, R2 (optional), p
  cols <- setdiff(names(df), "Level")
  n    <- nrow(df)
  xs   <- setNames(if (length(cols) > 2) c(0.44, 0.70, 1) else c(0.72, 1), cols)
  long <- do.call(rbind, lapply(seq_len(n), function(i)
    data.frame(x = c(0, xs), y = -i, label = c(as.character(df$Level[i]), unlist(df[i, cols])),
               h = c(0, rep(1, length(cols))), stringsAsFactors = FALSE)))
  head <- data.frame(x = c(0, xs), y = 0.2, label = c("Level", cols), h = c(0, rep(1, length(cols))))
  hl_y <- -which(as.character(df$Level) == focus_level)
  ggplot() +
    annotation_custom(grid::roundrectGrob(r = grid::unit(0.3, "snpc"), gp = grid::gpar(fill = hl_fill, col = hl_col, lwd = 1.2)),
                      xmin = -0.05, xmax = 1.05, ymin = hl_y - 0.45, ymax = hl_y + 0.45) +
    geom_text(data = long, aes(x, y, label = label, hjust = h), size = 2.8) +
    geom_text(data = head, aes(x, y, label = label, hjust = h), size = 2.8, fontface = "bold") +
    annotate("segment", x = -0.05, xend = 1.05, y = -0.3, yend = -0.3, linewidth = 0.3) +
    scale_x_continuous(limits = c(-0.12, 1.12), expand = c(0, 0)) +
    scale_y_continuous(limits = c(-n - 0.8, 1.2), expand = c(0, 0)) +
    labs(title = title, tag = tag) + theme_void() + tag_theme +
    theme(plot.margin = margin(5.5, 5.5, 5.5, 5.5))
}

sr_tab <- suff %>% filter(index == "SR") %>%
  mutate(Level = factor(level, levels = rev(levels_tax))) %>% arrange(Level) %>%
  transmute(Level, r = sprintf("%.3f", r), `R²` = sprintf("%.3f", R2))
sh_tab <- suff %>% filter(index == "Shannon") %>%
  mutate(Level = factor(level, levels = rev(levels_tax))) %>% arrange(Level) %>%
  transmute(Level, r = sprintf("%.3f", r), `R²` = sprintf("%.3f", R2))
mt_tab <- mantel_tab %>% mutate(Level = factor(level, levels = rev(levels_tax))) %>% arrange(Level) %>%
  transmute(Level = as.character(Level), Stress = sprintf("%.3f", stress_all[as.character(Level)]), `R²` = sprintf("%.3f", r^2))

lev_nmds <- function(l) {
  d <- ord_all[[l]]
  ggplot(d, aes(NMDS1, NMDS2)) +
    geom_path(arrow = arrow(type = "closed", length = unit(0.02, "npc")), colour = "grey50") +
    geom_point(aes(colour = year), size = 1.4) +
    scale_colour_viridis_c(option = "mako", end = 0.9) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5, size = 2.6,
             label = if (l == "Species") sprintf("Stress = %.3f", stress_all[[l]]) else "") +
    labs(title = l, tag = next_tag()) +
    theme_paper + tag_theme + theme(legend.position = "none", axis.title = element_text(size = 7.5),
                                    axis.text = element_text(size = 7),
                                    plot.background = element_blank(),
                                    panel.background = element_rect(fill = if (l == focus_level) hl_fill else "white", colour = NA))
}

left  <- c(lapply(c("Class", "Order", "Family", "Genus"), suff_panel),
           list(res_table(sr_tab, "Species richness", next_tag()), res_table(sh_tab, "Shannon index", next_tag())))
right <- c(lapply(c("Class", "Order", "Family", "Genus", "Species"), lev_nmds),
           list(res_table(mt_tab, "Mantel test (R²) and stress", next_tag())))
png(here("figures", "sufficiency.png"), width = 7.2, height = 5.6, units = "in", res = 300)
grid::grid.newpage()
grid::pushViewport(grid::viewport(layout = grid::grid.layout(3, 4)))
place <- function(p, i, col0, tag, hl = FALSE)
  draw_cell(p, tag, (i - 1) %/% 2 + 1, col0 + (i - 1) %% 2, hl = hl, plot_h = 0.93, card_h = 0.93, card_w = 0.94)
for (i in seq_along(left)) place(left[[i]], i, 1, letters[i], hl = i == 3)           # Family scatter (c)
for (i in seq_along(right)) place(right[[i]], i, 3, letters[6 + i], hl = i == 3)     # Family nMDS (i)
grid::grid.lines(x = grid::unit(0.5, "npc"), y = grid::unit(c(0.02, 0.98), "npc"),
                 gp = grid::gpar(col = "grey55", lwd = 0.8))      # separator between the two grids
invisible(dev.off())

# ---- Surrogacy: taxonomic groups as subsets ------------------------------------
groups <- list(Polychaeta = BDD$Class == "Polychaeta", Mollusca = BDD$Phylum == "Mollusca",
               Amphipoda = BDD$Order == "Amphipoda")
surr <- list()
for (g in names(groups)) {
  sel <- groups[[g]]; Ms <- S_full[sel, , drop = FALSE]
  sub <- idx_of(Ms)
  for (ix in c("SR", "Shannon")) {
    fit <- lm(species_idx[[ix]] ~ sub[[ix]]); sm <- summary(fit)
    surr[[length(surr) + 1]] <- data.frame(group = g, species = sum(sel), index = ix,
      r = round(sqrt(sm$r.squared) * sign(coef(fit)[2]), 3), R2 = round(sm$r.squared, 3),
      p = signif(coef(sm)[2, 4], 3))
  }
  mt <- mantel(bray, vegdist(t(log1p(Ms)), method = "bray"), method = "pearson", permutations = 999)
  surr[[length(surr) + 1]] <- data.frame(group = g, species = sum(sel), index = "Bray-Curtis (Mantel)",
    r = round(mt$statistic, 3), R2 = round(mt$statistic^2, 3), p = signif(mt$signif, 3))
}
surr <- do.call(rbind, surr); print(surr)

# Surrogacy figure: Amphipoda (best surrogate) against the whole community for SR, Shannon and Bray-Curtis
surr_tag  <- 0
next_tag  <- function() { surr_tag <<- surr_tag + 1; letters[surr_tag] }

scatter_panel <- function(df, xlab, ylab, title) {
  fit <- cor.test(df$x, df$y, method = "pearson")
  lab <- sprintf("r = %.2f, R\u00b2 = %.2f", unname(fit$estimate), unname(fit$estimate)^2)
  ggplot(df, aes(x, y)) +
    geom_point(aes(colour = year), size = 1.5) +
    geom_smooth(method = "lm", se = TRUE, level = 0.95, colour = "grey40", fill = "grey60", alpha = 0.25, linewidth = 0.5, formula = y ~ x) +
    annotate("text", x = -Inf, y = Inf, label = lab, hjust = -0.05, vjust = 1.4, size = 2.5) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
    scale_colour_viridis_c(option = "mako", end = 0.9) +
    labs(x = xlab, y = ylab, title = title, tag = next_tag()) + theme_paper +
    theme(legend.position = "none", aspect.ratio = 1, plot.tag = element_text(face = "bold", size = 11),
          plot.title = element_text(face = "bold", size = 8.5, hjust = 0.5),
          axis.title = element_text(size = 7.5), axis.text = element_text(size = 7))
}
g      <- "Amphipoda"                                                  # scatter plots for the best surrogate only
am     <- S_full[groups[[g]], , drop = FALSE]
am_idx <- idx_of(am)
nmds_sq <- function(scores_df, stress, title) {                              # square nMDS panel coloured by year, stress on plot
  d <- data.frame(NMDS1 = scores_df$NMDS1, NMDS2 = scores_df$NMDS2, year = years)
  ggplot(d, aes(NMDS1, NMDS2)) +
    geom_path(arrow = arrow(type = "closed", length = unit(0.02, "npc")), colour = "grey50") +
    geom_point(aes(colour = year), size = 1.5) +
    scale_colour_viridis_c(option = "mako", end = 0.9) +
    annotate("text", x = -Inf, y = Inf, label = sprintf("Stress = %.3f", stress), hjust = -0.05, vjust = 1.4, size = 2.5) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
    labs(x = "NMDS1", y = "NMDS2", title = title, tag = next_tag()) + theme_paper +
    theme(legend.position = "none", aspect.ratio = 1, plot.tag = element_text(face = "bold", size = 11),
          plot.title = element_text(face = "bold", size = 8.5, hjust = 0.5),
          axis.title = element_text(size = 7.5), axis.text = element_text(size = 6.5))
}
surr_panels <- list(
  scatter_panel(data.frame(x = species_idx$SR, y = am_idx$SR, year = years), "Whole-community SR", paste(g, "SR"), title = "Species richness"),
  scatter_panel(data.frame(x = species_idx$Shannon, y = am_idx$Shannon, year = years), "Whole-community Shannon", paste(g, "Shannon"), title = "Shannon index"),
  nmds_sq(sc[, c("NMDS1", "NMDS2")], nmds_bray$stress, "All species"),        # whole community, same ordination as the beta figure
  {                                                                           # Amphipoda species matrix (Bray-Curtis, log(x + 1))
    d_am <- vegdist(t(log1p(am)), method = "bray"); set.seed(1)
    nm_am <- metaMDS(d_am, k = 2, autotransform = FALSE, trace = FALSE)
    cat("Amphipoda nMDS stress:", round(nm_am$stress, 3), "\n")
    nmds_sq(as.data.frame(scores(nm_am, display = "sites")), nm_am$stress, "Amphipoda")
  })
png(here("figures", "surrogacy.png"), width = 7.2, height = 2.2, units = "in", res = 300)
grid::grid.newpage()
grid::pushViewport(grid::viewport(layout = grid::grid.layout(1, 4)))
for (i in seq_along(surr_panels))
  print(surr_panels[[i]], vp = grid::viewport(layout.pos.row = 1, layout.pos.col = i))
invisible(dev.off())

cat("Done. Figures written to figures/.\n")
