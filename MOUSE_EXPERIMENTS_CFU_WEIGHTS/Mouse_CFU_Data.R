############################################################
# Mouse CFU + Weight Analysis
# - Feces CFU: Wilcoxon per day; trend (replicate-level LM)
# - Organ CFU: Wilcoxon per organ
# - Weights & Gains: means ± SD, Wilcoxon per day
############################################################

# =========================
# Section 0. Packages
# =========================
ensure_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

ensure_pkg("dplyr")
ensure_pkg("tibble")
ensure_pkg("tidyr")
ensure_pkg("broom")

# Optional exports (uncomment if needed)
# ensure_pkg("writexl")
# ensure_pkg("officer")
# ensure_pkg("flextable")

cat("\nPackage versions:\n")
cat(sprintf("R:        %s\n", R.version.string))
for (p in c("dplyr","tibble","tidyr","broom")) {
  cat(sprintf("%-8s %s\n", paste0(p, ":"), as.character(packageVersion(p))))
}
cat("\n")

# =========================
# Section 1. Data (raw replicates)
# =========================
# Feces CFU (Control & Lysin), Days 21–24, n=6 per group/day
feces_control <- list(
  D21 = c(402.7, 937.2, 2040.2, 2125.6, 3599.9, 1366.9),
  D22 = c(1131.3, 334.9, 5154.7, 2544.1, 1818.7, 1203.8),
  D23 = c(615.3, 167.3, 1963.0, 1065.2, 4200.6, 1831.6),
  D24 = c(571.4, 514.3, 908.6, 1039.7, 3073.2, 2693.4)
)
feces_lysin <- list(
  D21 = c(1027.7, 1574.1, 859.4, 4091.6, 259.2, 1021.2),
  D22 = c(653.3, 965.1, 756.9, 3915.6, 1907.9, 802.8),
  D23 = c(258.2, 436.8, 634.0, 713.4, 160.7, 530.5),
  D24 = c(602.2, 482.5, 889.0, 742.8, 254.0, 647.8)
)

# Organ CFU (Control & Lysin), 6 replicates per group
organs_control <- list(
  Duodenum = c(295.9, 127.3, 197.6, 177.4, 267.5, 181.5),
  Jejunum  = c(230.3, 422.8, 270.0, 236.6, 673.8, 2214.7),
  Ileum    = c(41.1, 258.3, 2399.1, 796.3, 608.9, 669.0)
)
organs_lysin <- list(
  Duodenum = c(75.8, 154.4, 55.6, 287.2, 176.6, 468.3),
  Jejunum  = c(597.0, 176.6, 241.9, 211.9, 155.8, 80.8),
  Ileum    = c(242.4, 194.2, 231.3, 297.5, 96.1, 362.7)
)

# Mouse weights (IDs preserved; decimal points must be '.')
weights <- tibble::tribble(
  ~MouseID, ~Treatment,    ~Day0, ~Day7, ~Day14, ~Day21, ~Day22, ~Day23, ~Day24, ~Gain21, ~Gain22, ~Gain23, ~Gain24,
  "XB_0119", "+ Lysin",     10.2,  11.6,  12.6,   14.0,   14.6,   14.7,   14.7,    3.8,     4.4,     4.5,     4.5,
  "XB_0120", "+ Lysin",     10.4,  11.9,  13.4,   15.2,   15.5,   15.9,   16.0,    4.8,     5.1,     5.5,     5.6,
  "XB_0121", "PBS control", 10.4,  12.0,  13.0,   14.3,   14.8,   14.9,   14.8,    3.9,     4.4,     4.5,     4.4,
  "XB_0122", "PBS control",  9.6,  11.2,  12.2,   14.3,   13.7,   13.5,   13.8,    4.7,     4.1,     3.9,     4.2,
  "XB_0123", "PBS control", 10.2,  10.8,  12.2,   14.0,   14.0,   14.0,   14.3,    3.8,     3.8,     3.8,     4.1,
  "XB_0124", "PBS control",  9.9,  11.7,  12.3,   13.2,   14.1,   14.4,   14.3,    3.3,     4.2,     4.5,     4.4,
  "XB_0125", "+ Lysin",      9.6,  10.9,  12.5,   14.0,   14.3,   14.5,   14.0,    4.4,     4.7,     4.9,     4.4,
  "XB_0126", "+ Lysin",      9.7,  12.1,  13.2,   14.9,   15.2,   15.0,   14.9,    5.2,     5.5,     5.3,     5.2,
  "XB_0127", "+ Lysin",      9.7,  11.9,  13.1,   14.6,   15.0,   14.9,   15.1,    4.9,     5.3,     5.2,     5.4,
  "XB_0128", "+ Lysin",     10.4,  11.7,  13.6,   14.6,   15.1,   15.5,   15.6,    4.2,     4.7,     5.1,     5.2,
  "XB_0129", "PBS control",  9.9,  12.0,  13.9,   15.9,   16.1,   16.6,   17.1,    6.0,     6.2,     6.7,     7.2,
  "XB_0130", "PBS control",  9.6,  11.7,  13.8,   15.8,   16.3,   16.6,   16.9,    6.2,     6.7,     7.0,     7.3
)

# =========================
# Section 2.
# =========================
mean_sd <- function(x, digits = 2) paste0(round(mean(x), digits), " ± ", round(sd(x), digits))
fmt_p <- function(p) ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))

wilcox_two <- function(x, y) stats::wilcox.test(x, y, exact = FALSE)$p.value

# =========================
# Section 3. Feces CFU – Wilcoxon per day (Lysin vs Control)
# =========================
feces_tests <- lapply(names(feces_control), function(day) {
  x <- feces_control[[day]]
  y <- feces_lysin[[day]]
  p <- wilcox_two(x, y)
  tibble::tibble(
    Day = day,
    Control_mean_SD = mean_sd(x),
    Lysin_mean_SD   = mean_sd(y),
    p_value         = p,
    `p-value`       = fmt_p(p),
    Test            = "Wilcoxon rank-sum (two-sided)"
  )
}) |> dplyr::bind_rows()

cat("=== Feces CFU (Control vs Lysin) by day ===\n")
print(feces_tests)

# =========================
# Section 4. Organ CFU – Wilcoxon per organ
# =========================
organ_tests <- lapply(names(organs_control), function(org) {
  x <- organs_control[[org]]
  y <- organs_lysin[[org]]
  p <- wilcox_two(x, y)
  tibble::tibble(
    Organ           = org,
    Control_mean_SD = mean_sd(x),
    Lysin_mean_SD   = mean_sd(y),
    N_per_group     = length(x),
    p_value         = p,
    `p-value`       = fmt_p(p),
    Test            = "Wilcoxon rank-sum (two-sided)"
  )
}) |> dplyr::bind_rows()

cat("\n=== Organ CFU (Control vs Lysin) ===\n")
print(organ_tests)

# =========================
# Section 5. Feces CFU – Trend (replicate-level LM)
# (We avoid funtimes::notrend_test due to only 4 time points)
# =========================
feces_all <- tibble::tibble(
  Day   = rep(21:24, each = 6, times = 2),
  Group = rep(c("Control", "Lysin"), each = 24),
  CFU   = c(unlist(feces_control), unlist(feces_lysin))
)

lm_ctrl <- stats::lm(CFU ~ Day, data = dplyr::filter(feces_all, Group == "Control"))
lm_lys  <- stats::lm(CFU ~ Day, data = dplyr::filter(feces_all, Group == "Lysin"))

trend_ctrl <- broom::tidy(lm_ctrl) |> dplyr::filter(term == "Day") |>
  dplyr::transmute(Group = "Control",
                   Slope_CFU_per_day = round(estimate, 1),
                   SE = round(std.error, 1),
                   p_value = p.value,
                   `p-value` = fmt_p(p.value),
                   Model = "Linear regression (CFU ~ Day)")
trend_lys <- broom::tidy(lm_lys) |> dplyr::filter(term == "Day") |>
  dplyr::transmute(Group = "Lysin",
                   Slope_CFU_per_day = round(estimate, 1),
                   SE = round(std.error, 1),
                   p_value = p.value,
                   `p-value` = fmt_p(p.value),
                   Model = "Linear regression (CFU ~ Day)")

trend_results <- dplyr::bind_rows(trend_ctrl, trend_lys)

cat("\n=== Feces CFU Trend (replicate-level) ===\n")
print(trend_results)

# =========================
# Section 6. Mouse weights – group means ± SD and Wilcoxon per day
# =========================
# (A) Weights by day
weight_days <- c("Day0","Day7","Day14","Day21","Day22","Day23","Day24")
weights_summary <- lapply(weight_days, function(col) {
  ctrl <- dplyr::filter(weights, Treatment == "PBS control")[[col]]
  lys  <- dplyr::filter(weights, Treatment == "+ Lysin")[[col]]
  p <- wilcox_two(lys, ctrl)
  tibble::tibble(
    Timepoint       = col,
    Control_mean_SD = mean_sd(ctrl),
    Lysin_mean_SD   = mean_sd(lys),
    p_value         = p,
    `p-value`       = fmt_p(p),
    Test            = "Wilcoxon rank-sum (two-sided)"
  )
}) |> dplyr::bind_rows()

cat("\n=== Mouse weights (g): Control vs Lysin per day ===\n")
print(weights_summary)

# (B) Gains vs Day 0
gain_days <- c("Gain21","Gain22","Gain23","Gain24")
gains_summary <- lapply(gain_days, function(col) {
  ctrl <- dplyr::filter(weights, Treatment == "PBS control")[[col]]
  lys  <- dplyr::filter(weights, Treatment == "+ Lysin")[[col]]
  p <- wilcox_two(lys, ctrl)
  tibble::tibble(
    Timepoint       = col,
    Control_mean_SD = mean_sd(ctrl),
    Lysin_mean_SD   = mean_sd(lys),
    p_value         = p,
    `p-value`       = fmt_p(p),
    Test            = "Wilcoxon rank-sum (two-sided)"
  )
}) |> dplyr::bind_rows()

cat("\n=== Weight gains vs Day 0 (g): Control vs Lysin ===\n")
print(gains_summary)



############################################################
# Notes:
# - Wilcoxon tests are two-sided (Mann–Whitney).
# - Trend analysis uses replicate-level linear regression (4 timepoints).
# - Report exact p-values; consider multiple-testing correction if needed.
############################################################
