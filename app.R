# ============================================================
# Data Science for Chemists — Short Course Shiny App
# Chapters 1, 2, and 5
#
# Required packages (install once before running):
#   install.packages(c("shiny","bslib","tidyverse","DT",
#                      "mice","VIM","umap","cowplot"))
#
# Chapter 5: place NHANES.csv in this same folder.
# ============================================================

suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(tidyverse)
  library(DT)
  library(mice)
  library(VIM)
  library(umap)
  library(cowplot)
})

# ─────────────────────────────────────────────────────────────
# GLOBAL: Pre-compute all deterministic data at startup
# (set.seed(42) ensures every participant sees the same data)
# ─────────────────────────────────────────────────────────────

message("[1/6] Generating chemistry dataset...")
set.seed(42)
n_samples <- 60

samples <- tibble(
  SampleID   = sprintf("S%03d", 1:n_samples),
  Location   = sample(c("Upstream", "Midstream", "Downstream"), n_samples, TRUE),
  MatrixType = sample(c("Water", "Sediment"), n_samples, TRUE, prob = c(0.6, 0.4)),
  Depth_m    = ifelse(MatrixType == "Water",
                      round(runif(n_samples, 0.2, 2.0), 2),
                      round(runif(n_samples, 0.0, 0.3), 2)),
  Date       = as.Date("2026-01-01") + sample(0:55, n_samples, TRUE),
  TOC_level  = sample(c("Low", "Medium", "High"), n_samples, TRUE, prob = c(0.4, 0.4, 0.2)),
  QC_pass    = sample(c(TRUE, FALSE), n_samples, TRUE, prob = c(0.9, 0.1))
) %>%
  mutate(TOC_level = factor(TOC_level, levels = c("Low", "Medium", "High"), ordered = TRUE))

analytes       <- c("PFOS", "PFOA", "PFHxS", "PFNA", "PFUnDA", "PCB52", "PCB101", "Nap", "Phe")
toc_multiplier <- c(Low = 0.8, Medium = 1.0, High = 1.4)

conc_wide <- map_dfc(analytes, function(a) {
  base          <- runif(1, 0.05, 2.0)
  loc_effect    <- case_when(
    samples$Location == "Upstream"   ~ 0.8,
    samples$Location == "Midstream"  ~ 1.0,
    samples$Location == "Downstream" ~ 1.3
  )
  matrix_effect <- ifelse(samples$MatrixType == "Sediment", 1.8, 1.0)
  toc_effect    <- toc_multiplier[as.character(samples$TOC_level)]
  x             <- base * loc_effect * matrix_effect * toc_effect * exp(rnorm(n_samples, 0, 0.35))
  nd_prob       <- if (a %in% c("PFUnDA", "PCB101")) 0.25 else 0.15
  x[runif(n_samples) < nd_prob] <- NA_real_
  tibble(!!a := round(x, 3))
})

chem_wide <- bind_cols(samples, conc_wide)

RL     <- tibble(Analyte = analytes,
                 RL      = c(0.5, 0.4, 0.3, 0.2, 0.15, 0.05, 0.05, 0.02, 0.02))
RL_vec <- setNames(RL$RL, RL$Analyte)

chem_long <- chem_wide %>%
  pivot_longer(cols = all_of(analytes), names_to = "Analyte", values_to = "Concentration") %>%
  left_join(RL, by = "Analyte") %>%
  mutate(Detect   = !is.na(Concentration),
         Conc_imp = ifelse(is.na(Concentration), RL / 2, Concentration))

# RL/2 imputation (used as predictors in LM + Mahalanobis)
chem_rl2 <- chem_wide
for (a in analytes) chem_rl2[[a]] <- ifelse(is.na(chem_rl2[[a]]), RL_vec[a] / 2, chem_rl2[[a]])

# ── LM imputation ────────────────────────────────────────────
message("[2/6] Computing least-squares imputation...")
chem_lm   <- chem_wide
X_pred_rl <- chem_rl2 %>% select(all_of(analytes))

for (a in analytes) {
  y         <- chem_wide[[a]]
  miss_idx  <- which(is.na(y))
  if (length(miss_idx) == 0) next
  pred_vars <- setdiff(analytes, a)
  train_df  <- X_pred_rl %>% mutate(y = y) %>% filter(!is.na(y))
  test_df   <- X_pred_rl[miss_idx, pred_vars, drop = FALSE]
  fit       <- lm(y ~ ., data = train_df %>% select(all_of(pred_vars), y))
  chem_lm[[a]][miss_idx] <- pmax(predict(fit, newdata = test_df), 0)
}

# LM actual vs predicted (on observed rows, for diagnostics)
avp_lm <- map_dfr(analytes, function(a) {
  y         <- chem_wide[[a]]
  if (sum(!is.na(y)) < 5) return(NULL)
  pred_vars <- setdiff(analytes, a)
  train_df  <- X_pred_rl %>% mutate(y = y) %>% filter(!is.na(y))
  fit       <- lm(y ~ ., data = train_df %>% select(all_of(pred_vars), y))
  yhat      <- pmax(predict(fit, newdata = train_df %>% select(all_of(pred_vars))), 0)
  tibble(Analyte = a, Actual = train_df$y, Predicted = yhat, Residual = train_df$y - yhat)
})

# ── kNN imputation (full) ─────────────────────────────────────
message("[3/6] Computing kNN imputation (full)...")
knn_full   <- VIM::kNN(chem_wide %>% select(all_of(analytes)), k = 5, imp_var = FALSE)
chem_knn   <- chem_wide
chem_knn[analytes] <- knn_full[analytes]

# ── kNN holdout validation (10% mask) ─────────────────────────
message("[4/6] Computing kNN holdout validation...")
set.seed(42)
chem_masked <- chem_wide
holdout_idx <- map(analytes, function(a) {
  obs <- which(!is.na(chem_wide[[a]]))
  sample(obs, max(1, floor(length(obs) * 0.10)))
})
names(holdout_idx) <- analytes

truth_long <- map_dfr(analytes, function(a) {
  idx <- holdout_idx[[a]]
  tibble(SampleID = chem_wide$SampleID[idx], Analyte = a, Actual = chem_wide[[a]][idx])
})
for (a in analytes) chem_masked[[a]][holdout_idx[[a]]] <- NA_real_

knn_masked <- VIM::kNN(chem_masked %>% select(all_of(analytes)), k = 5, imp_var = FALSE)
chem_knn_masked           <- chem_masked
chem_knn_masked[analytes] <- knn_masked[analytes]

pred_long_knn <- chem_knn_masked %>%
  select(SampleID, all_of(analytes)) %>%
  pivot_longer(all_of(analytes), names_to = "Analyte", values_to = "Predicted")

avp_knn <- truth_long %>%
  left_join(pred_long_knn, by = c("SampleID", "Analyte")) %>%
  mutate(Residual = Actual - Predicted)

# ── MICE imputation ───────────────────────────────────────────
message("[5/6] Computing MICE imputation (~15 s)...")
imp_mice        <- mice(chem_wide %>% select(all_of(analytes)),
                        m = 5, method = "pmm", maxit = 10, seed = 42, printFlag = FALSE)
chem_mice       <- chem_wide
chem_mice[analytes] <- complete(imp_mice, action = 1)

# ── Mahalanobis distances ─────────────────────────────────────
message("[6/6] Computing Mahalanobis distances...")
X_rl2_mat  <- chem_rl2 %>% select(all_of(analytes)) %>% as.matrix()
md_dist    <- mahalanobis(X_rl2_mat, colMeans(X_rl2_mat), cov(X_rl2_mat))
p_md       <- pchisq(md_dist, df = ncol(X_rl2_mat), lower.tail = FALSE)

md_results <- chem_wide %>%
  select(SampleID, Location, MatrixType) %>%
  mutate(Mahalanobis = round(md_dist, 2),
         p_value     = round(p_md, 4),
         MV_outlier  = p_md < 0.001) %>%
  arrange(p_value)

message("App ready!")

# ── Load NHANES (Chapter 5) ───────────────────────────────────
nhanes <- tryCatch(read.csv("NHANES.csv"), error = function(e) NULL)

# ─────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────
iqr_flags <- function(x) {
  q1  <- quantile(x, 0.25, na.rm = TRUE)
  q3  <- quantile(x, 0.75, na.rm = TRUE)
  iqr <- q3 - q1
  out <- (x < q1 - 1.5 * iqr) | (x > q3 + 1.5 * iqr)
  out[is.na(x)] <- NA
  out
}

code_box <- function(txt) {
  tags$pre(
    style = paste("background:#f8f9fa; border:1px solid #dee2e6; border-radius:6px;",
                  "padding:14px; font-size:12.5px; line-height:1.55; overflow-x:auto;"),
    tags$code(class = "language-r", txt)
  )
}

dist_plot <- function(x, title, fill_col) {
  sw <- shapiro.test(x)
  tibble(v = x) %>%
    ggplot(aes(v)) +
    geom_histogram(bins = 20, fill = fill_col, color = "white") +
    labs(title = paste0(title, "\nShapiro p = ", round(sw$p.value, 3)),
         x = NULL, y = "Count") +
    theme_bw(base_size = 12)
}

# ─────────────────────────────────────────────────────────────
# UI
# ─────────────────────────────────────────────────────────────
ui <- navbarPage(
  title = "Data Science for Chemists — Short Course",
  theme = bs_theme(bootswatch = "flatly"),

  # ══════════════════════════════════════════════════════════
  # CHAPTER 1: Introduction
  # ══════════════════════════════════════════════════════════
  tabPanel(
    "Ch 1: Introduction",
    fluidPage(
      tabsetPanel(type = "tabs",

        # ── 1. Dataset Generation ──────────────────────────
        tabPanel("1. Dataset Generation", br(),
          fluidRow(
            column(5,
              h4("Generating a Realistic Chemistry Dataset"),
              p("60 samples · 9 analytes (PFAS/PCBs/PAHs) · structured metadata · ",
                "simulated non-detects driven by matrix type, location, and TOC level."),
              code_box(
'set.seed(42)
n_samples <- 60

samples <- tibble(
  SampleID   = sprintf("S%03d", 1:n_samples),
  Location   = sample(c("Upstream","Midstream","Downstream"),
                      n_samples, replace = TRUE),
  MatrixType = sample(c("Water","Sediment"), n_samples,
                      TRUE, prob = c(0.6, 0.4)),
  Depth_m    = ifelse(MatrixType == "Water",
                 round(runif(n_samples, 0.2, 2.0), 2),
                 round(runif(n_samples, 0.0, 0.3), 2)),
  Date       = as.Date("2026-01-01") +
                 sample(0:55, n_samples, TRUE),
  TOC_level  = factor(
                 sample(c("Low","Medium","High"), n_samples,
                        TRUE, prob = c(0.4, 0.4, 0.2)),
                 levels = c("Low","Medium","High"),
                 ordered = TRUE),
  QC_pass    = sample(c(TRUE,FALSE), n_samples,
                      TRUE, prob = c(0.9, 0.1))
)

analytes <- c("PFOS","PFOA","PFHxS","PFNA","PFUnDA",
              "PCB52","PCB101","Nap","Phe")

# Concentrations: base * location * matrix * TOC effects
#   + log-normal noise
# Non-detects: 15-25% probability per analyte -> NA'
              )
            ),
            column(7,
              h5("chem_wide (first 20 rows):"),
              DTOutput("ch1_tbl"),
              br(),
              h5("glimpse(chem_wide):"),
              verbatimTextOutput("ch1_glimpse")
            )
          )
        ),

        # ── 2. Data Structures ─────────────────────────────
        tabPanel("2. Data Structures", br(),
          h4("Vectors, Matrices, Arrays, Data Frames"),
          fluidRow(
            column(6,
              h5("Vector (1D)"),
              code_box(
'pfos_vec <- chem_wide$PFOS
length(pfos_vec)
head(pfos_vec)'
              ),
              verbatimTextOutput("ch1_vector"),
              br(),
              h5("Data Frame / Tibble (2D, mixed types)"),
              code_box(
'# Holds numeric, character, factor,
# date, logical in one object
is.data.frame(chem_wide)
map_chr(chem_wide, ~ class(.x)[1])'
              ),
              verbatimTextOutput("ch1_df")
            ),
            column(6,
              h5("Matrix (2D, numeric only)"),
              code_box(
'X <- chem_wide %>%
  select(all_of(analytes)) %>%
  as.matrix()
dim(X)   # rows = samples, cols = analytes
X[1:3, 1:4]'
              ),
              verbatimTextOutput("ch1_matrix"),
              br(),
              h5("Array (3D: Sample × Analyte × Version)"),
              code_box(
'A <- array(NA_real_,
  dim      = c(60, 9, 2),
  dimnames = list(
    chem_wide$SampleID,
    analytes,
    c("Observed", "Imputed")
  ))
dim(A)
A[1:3, 1:3, ]   # 3-way peek'
              ),
              verbatimTextOutput("ch1_array")
            )
          )
        ),

        # ── 3. Matrix Operations ───────────────────────────
        tabPanel("3. Matrix Operations", br(),
          h4("Addition, Subtraction, Multiplication"),
          fluidRow(
            column(5,
              code_box(
'# Use RL/2-imputed matrix (no NAs)
X_imp <- A[,, "Imputed"]

# Addition (commutative: A+B == B+A)
noise <- matrix(rnorm(60*9, 0, 0.01), nrow = 60)
all.equal(X_imp + noise, noise + X_imp)

# Subtraction (NOT commutative)
diff1 <- X_imp - noise
diff2 <- noise - X_imp
mean(abs(diff1 - diff2))  # > 0

# Matrix multiply:
#   (60x9) %*% (9x1) -> (60x1)
# Weighted contamination score
weights <- runif(9, 0.5, 1.5)
weights <- weights / sum(weights)
W       <- matrix(weights, ncol = 1)
score   <- X_imp %*% W
dim(score)   # 60 x 1
head(score)'
              )
            ),
            column(7,
              verbatimTextOutput("ch1_matops")
            )
          )
        ),

        # ── 4. Data Manipulation ───────────────────────────
        tabPanel("4. Data Manipulation", br(),
          h4("Subset, Split, Pivot (Wide ↔ Long), Join"),
          fluidRow(
            column(5,
              h5("Subset vs Split"),
              code_box(
'# Subset: filter rows
sediment_only <- chem_wide %>%
  filter(MatrixType == "Sediment")
nrow(sediment_only)

# Split: divide into a list by group
by_location <- split(chem_wide,
                     chem_wide$Location)
lapply(by_location, nrow)'
              ),
              verbatimTextOutput("ch1_subset"),
              br(),
              h5("Joins"),
              code_box(
'# left_join: keep ALL left rows
chem_joined <- chem_wide %>%
  left_join(field, by = "SampleID")

# inner_join: only matching rows
# full_join:  all rows from both tables
# anti_join:  left rows NOT in right

# Update QC flags from a new table
chem_updated <- chem_wide %>%
  left_join(qc_updates, by = "SampleID") %>%
  mutate(
    QC_pass = ifelse(!is.na(QC_pass_new),
                     QC_pass_new, QC_pass)
  ) %>%
  select(-QC_pass_new)'
              ),
              verbatimTextOutput("ch1_join")
            ),
            column(7,
              h5("Wide → Long (tidy format for analysis)"),
              code_box(
'chem_long <- chem_wide %>%
  pivot_longer(
    cols      = all_of(analytes),
    names_to  = "Analyte",
    values_to = "Concentration"
  ) %>%
  left_join(RL, by = "Analyte") %>%
  mutate(
    Detect   = !is.na(Concentration),
    Conc_imp = ifelse(is.na(Concentration),
                      RL / 2, Concentration)
  )

# Detection frequency by analyte
chem_long %>%
  group_by(Analyte) %>%
  summarise(
    DetectRate = mean(Detect),
    Median_Imp = median(Conc_imp)
  ) %>%
  arrange(desc(DetectRate))'
              ),
              verbatimTextOutput("ch1_pivot")
            )
          )
        ),

        # ── 5. Plots ───────────────────────────────────────
        tabPanel("5. Plots", br(),
          h4("Visualizing Chemical Data with ggplot2"),
          fluidRow(
            column(6,
              h5("Detection Rate by Analyte"),
              plotOutput("ch1_detect_plot", height = "320px"),
              br(),
              code_box(
'chem_long %>%
  group_by(Analyte) %>%
  summarise(DetectRate = mean(Detect)) %>%
  ggplot(aes(x = reorder(Analyte, DetectRate),
             y = DetectRate)) +
  geom_col() +
  coord_flip() +
  labs(title = "Detection Rate by Analyte",
       x = NULL, y = "Detect rate")'
              )
            ),
            column(6,
              h5("Concentrations by Matrix Type"),
              plotOutput("ch1_box_plot", height = "320px"),
              br(),
              code_box(
'chem_long %>%
  ggplot(aes(x = MatrixType, y = Conc_imp)) +
  geom_boxplot(outlier.alpha = 0.3) +
  facet_wrap(~ Analyte, scales = "free_y") +
  labs(title = "Imputed Concentrations by Matrix Type",
       x = NULL, y = "Concentration (imputed)")'
              )
            )
          )
        )
      )
    )
  ),

  # ══════════════════════════════════════════════════════════
  # CHAPTER 2: Pre-processing
  # ══════════════════════════════════════════════════════════
  tabPanel(
    "Ch 2: Pre-processing",
    fluidPage(
      tabsetPanel(type = "tabs",

        # ── 1. IQR Outliers ────────────────────────────────
        tabPanel("1. Outliers — IQR", br(),
          fluidRow(
            column(4,
              h4("Univariate Outliers (IQR Method)"),
              selectInput("iqr_analyte", "Select analyte:",
                          choices = analytes, selected = "PFOS"),
              code_box(
'iqr_outlier_flags <- function(x) {
  q1  <- quantile(x, 0.25, na.rm = TRUE)
  q3  <- quantile(x, 0.75, na.rm = TRUE)
  iqr <- q3 - q1
  out <- (x < q1 - 1.5 * iqr) |
         (x > q3 + 1.5 * iqr)
  out[is.na(x)] <- NA
  out
}'
              ),
              br(),
              h5("Outlier count by analyte:"),
              DTOutput("iqr_summary")
            ),
            column(8,
              h5("Boxplot with IQR outliers highlighted (red ▲):"),
              plotOutput("iqr_plot", height = "430px")
            )
          )
        ),

        # ── 2. Mahalanobis ─────────────────────────────────
        tabPanel("2. Outliers — Mahalanobis", br(),
          fluidRow(
            column(5,
              h4("Multivariate Outliers (Mahalanobis Distance)"),
              code_box(
'# RL/2-filled matrix (no NAs required)
X_imp  <- chem_rl2 %>%
  select(all_of(analytes)) %>%
  as.matrix()

center <- colMeans(X_imp)
covmat <- cov(X_imp)

md   <- mahalanobis(X_imp, center, covmat)
p_md <- pchisq(md, df = ncol(X_imp),
               lower.tail = FALSE)
# p < 0.001 -> multivariate outlier'
              ),
              br(),
              h5("Top samples by Mahalanobis distance:"),
              DTOutput("mahal_tbl")
            ),
            column(7,
              plotOutput("mahal_bar", height = "380px"),
              plotOutput("fingerprint_plot", height = "280px")
            )
          )
        ),

        # ── 3. Missing Values ──────────────────────────────
        tabPanel("3. Missing Values", br(),
          fluidRow(
            column(5,
              h4("Exploring NA Patterns"),
              h5("By analyte:"),
              DTOutput("na_analyte_tbl"),
              br(),
              h5("By sample (top 10):"),
              DTOutput("na_sample_tbl")
            ),
            column(7,
              h5("Missingness map (red = non-detect / NA):"),
              plotOutput("na_heatmap", height = "460px"),
              br(),
              code_box(
'chem_wide %>%
  select(SampleID, all_of(analytes)) %>%
  pivot_longer(-SampleID,
    names_to  = "Analyte",
    values_to = "Conc") %>%
  mutate(Missing = is.na(Conc)) %>%
  ggplot(aes(x = Analyte, y = SampleID,
             fill = Missing)) +
  geom_tile()'
              )
            )
          )
        ),

        # ── 4. Imputation ──────────────────────────────────
        tabPanel("4. Imputation", br(),
          fluidRow(
            column(3,
              h4("Four Imputation Methods"),
              selectInput("imp_analyte", "Analyte for comparison:",
                          choices = analytes, selected = "PFOS"),
              selectInput("imp_method", "Show diagnostics for:",
                          choices = c("Least-squares (LM)"        = "lm",
                                      "kNN (holdout validation)"  = "knn"),
                          selected = "lm"),
              br(),
              tags$ul(
                tags$li(tags$b("RL/2"), " — replace ND with ½ reporting limit"),
                tags$li(tags$b("LM"), " — linear regression on other analytes"),
                tags$li(tags$b("kNN"), " — k-nearest neighbors (k = 5)"),
                tags$li(tags$b("MICE"), " — multiple imputation (PMM)")
              ),
              br(),
              h5("Imputed values for missing rows:"),
              DTOutput("imp_compare_tbl")
            ),
            column(9,
              h5("Actual vs Predicted (model fit diagnostic):"),
              plotOutput("imp_avp_plot", height = "380px"),
              br(),
              h5("Residuals (Actual − Predicted):"),
              plotOutput("imp_resid_plot", height = "300px")
            )
          )
        ),

        # ── 5. Feature Scaling ─────────────────────────────
        tabPanel("5. Feature Scaling", br(),
          fluidRow(
            column(3,
              h4("Scaling & Normality"),
              selectInput("scale_analyte", "Compound:",
                          choices = analytes, selected = "PFOS"),
              code_box(
'# Transformations
log10(x + eps)
(x - min(x)) / (max(x) - min(x))  # min-max
(x / sum(x)) * 100                # percent

# Shapiro-Wilk normality test
shapiro.test(x)
# H0: data is normal
# p > 0.05 -> fail to reject normality
# "Best" = highest p-value'
              ),
              br(),
              h5("Shapiro-Wilk results (sorted by p):"),
              DTOutput("shapiro_tbl")
            ),
            column(9,
              fluidRow(
                column(6,
                  h5("Raw"),
                  plotOutput("dist_raw", height = "210px")),
                column(6,
                  h5("log₁₀"),
                  plotOutput("dist_log", height = "210px"))
              ),
              fluidRow(
                column(6,
                  h5("Min-Max [0, 1]"),
                  plotOutput("dist_minmax", height = "210px")),
                column(6,
                  h5("Percent (%)"),
                  plotOutput("dist_percent", height = "210px"))
              )
            )
          )
        )
      )
    )
  ),

  # ══════════════════════════════════════════════════════════
  # CHAPTER 5: PCA vs UMAP
  # ══════════════════════════════════════════════════════════
  tabPanel(
    "Ch 5: PCA vs UMAP",
    fluidPage(
      fluidRow(
        column(3,
          wellPanel(
            h4("UMAP Parameters"),
            sliderInput("umap_k",    "n_neighbors:", min = 5,   max = 50,  value = 15,  step = 1),
            sliderInput("umap_dist", "min_dist:",    min = 0.0, max = 1.0, value = 0.3, step = 0.05),
            br(),
            actionButton("run_umap", "Run PCA + UMAP",
                         class = "btn-primary", width = "100%"),
            br(), br(),
            code_box(
'pca_res <- prcomp(X_raw,
  center = FALSE, scale. = FALSE)

cfg <- umap.defaults
cfg$n_neighbors  <- 15
cfg$min_dist     <- 0.3
cfg$n_components <- 2
umap_res <- umap(X_raw, config = cfg)

# Color by demographic variable
ggplot(umap_df,
  aes(UMAP1, UMAP2, color = Age)) +
  geom_point() +
  scale_color_brewer(palette = "Set2")'
            ),
            br(),
            p(class = "text-muted small",
              "Requires ", tags$code("NHANES.csv"), " in the app folder.",
              " Analytes: columns 7–41. Color variable: ", tags$code("Age"), " column.")
          )
        ),
        column(9,
          uiOutput("ch5_content")
        )
      )
    )
  )
)

# ─────────────────────────────────────────────────────────────
# SERVER
# ─────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  # ══════════════════════════════════════════════════════════
  # Chapter 1
  # ══════════════════════════════════════════════════════════

  output$ch1_tbl <- renderDT({
    chem_wide %>%
      mutate(Date = as.character(Date)) %>%
      head(20) %>%
      datatable(options = list(pageLength = 8, scrollX = TRUE), rownames = FALSE)
  })

  output$ch1_glimpse <- renderPrint({ glimpse(chem_wide) })

  output$ch1_vector <- renderPrint({
    v <- chem_wide$PFOS
    cat("length:", length(v), "\n")
    cat("head:  ", head(v, 8), "\n")
  })

  output$ch1_matrix <- renderPrint({
    X <- chem_wide %>% select(all_of(analytes)) %>% as.matrix()
    cat("dim(X):", dim(X), "  (rows = samples, cols = analytes)\n\n")
    cat("X[1:3, 1:4]:\n")
    print(round(X[1:3, 1:4], 3))
  })

  output$ch1_array <- renderPrint({
    A <- array(NA_real_,
               dim      = c(n_samples, length(analytes), 2),
               dimnames = list(chem_wide$SampleID, analytes, c("Observed", "Imputed")))
    for (j in seq_along(analytes)) {
      a   <- analytes[j]
      obs <- chem_wide[[a]]
      imp <- obs; imp[is.na(imp)] <- RL_vec[a] / 2
      A[, j, 1] <- obs; A[, j, 2] <- imp
    }
    cat("dim(A):", dim(A), "  (samples × analytes × version)\n\n")
    cat("A[1:3, 1:3, ]:\n")
    print(round(A[1:3, 1:3, ], 3))
  })

  output$ch1_df <- renderPrint({
    cat("class:", paste(class(chem_wide), collapse = ", "), "\n\n")
    cat("Column types:\n")
    print(map_chr(chem_wide, ~ class(.x)[1]))
  })

  output$ch1_matops <- renderPrint({
    set.seed(42)
    X_imp <- chem_rl2 %>% select(all_of(analytes)) %>% as.matrix()
    noise <- matrix(rnorm(n_samples * length(analytes), 0, 0.01), nrow = n_samples)

    cat("Addition commutative?",
        all.equal(X_imp + noise, noise + X_imp), "\n\n")

    d1 <- X_imp - noise; d2 <- noise - X_imp
    cat("Subtraction NOT commutative.\n")
    cat("  mean |diff1 - diff2|:", round(mean(abs(d1 - d2)), 5), "\n\n")

    weights <- runif(length(analytes), 0.5, 1.5)
    weights <- weights / sum(weights)
    score   <- X_imp %*% matrix(weights, ncol = 1)
    cat("score = X_imp %*% W\n")
    cat("  dim(score):", dim(score), "  (60 samples × 1)\n")
    cat("  head(score):", round(head(score, 5), 4), "\n")
  })

  output$ch1_subset <- renderPrint({
    sed    <- chem_wide %>% filter(MatrixType == "Sediment")
    by_loc <- split(chem_wide, chem_wide$Location)
    cat("Sediment rows:", nrow(sed), "of", nrow(chem_wide), "\n\n")
    cat("Rows by location:\n")
    print(sapply(by_loc, nrow))
  })

  output$ch1_pivot <- renderPrint({
    dr <- chem_long %>%
      group_by(Analyte) %>%
      summarise(DetectRate = round(mean(Detect), 3),
                Median_Imp = round(median(Conc_imp), 3),
                .groups    = "drop") %>%
      arrange(desc(DetectRate))
    print(as.data.frame(dr))
  })

  output$ch1_join <- renderPrint({
    set.seed(99)
    field <- chem_wide %>%
      select(SampleID, MatrixType) %>%
      mutate(pH         = round(rnorm(n(), mean = ifelse(MatrixType == "Water", 7.4, 7.0), sd = 0.3), 2),
             Conduct_uS = round(rlnorm(n(), meanlog = log(400), sdlog = 0.35), 1),
             Temp_C     = round(rnorm(n(), mean = 9, sd = 3), 1)) %>%
      select(SampleID, pH, Conduct_uS, Temp_C)

    field_missing <- field %>% slice(-c(1, 2))
    cat("inner_join rows (2 samples removed):", nrow(chem_wide %>% inner_join(field_missing, by = "SampleID")), "\n")
    cat("left_join rows  (all samples kept): ", nrow(chem_wide %>% left_join(field_missing,  by = "SampleID")), "\n")
    cat("\nSamples S001, S002 after left_join (pH = NA):\n")
    print(chem_wide %>% left_join(field_missing, by = "SampleID") %>%
            select(SampleID, pH) %>% head(3))
  })

  output$ch1_detect_plot <- renderPlot({
    chem_long %>%
      group_by(Analyte) %>%
      summarise(DetectRate = mean(Detect), .groups = "drop") %>%
      ggplot(aes(x = reorder(Analyte, DetectRate), y = DetectRate)) +
      geom_col(fill = "#2c7bb6") +
      coord_flip() +
      scale_y_continuous(limits = c(0, 1)) +
      labs(title = "Detection Rate by Analyte", x = NULL, y = "Detect rate") +
      theme_bw(base_size = 14)
  })

  output$ch1_box_plot <- renderPlot({
    chem_long %>%
      ggplot(aes(x = MatrixType, y = Conc_imp, fill = MatrixType)) +
      geom_boxplot(outlier.alpha = 0.3) +
      facet_wrap(~ Analyte, scales = "free_y") +
      labs(title = "Imputed Concentrations by Matrix Type",
           x = NULL, y = "Concentration (RL/2 imputed)") +
      theme_bw(base_size = 11) +
      theme(legend.position = "none")
  })

  # ══════════════════════════════════════════════════════════
  # Chapter 2
  # ══════════════════════════════════════════════════════════

  output$iqr_summary <- renderDT({
    tibble(
      Analyte   = analytes,
      N_outlier = map_int(analytes, ~ sum(iqr_flags(chem_wide[[.x]]), na.rm = TRUE)),
      N_nonNA   = map_int(analytes, ~ sum(!is.na(chem_wide[[.x]])))
    ) %>%
      mutate(Rate = round(N_outlier / N_nonNA, 3)) %>%
      arrange(desc(Rate)) %>%
      datatable(options = list(dom = "t", pageLength = 9), rownames = FALSE)
  })

  output$iqr_plot <- renderPlot({
    a <- input$iqr_analyte
    chem_wide %>%
      transmute(MatrixType, Location,
                Concentration = .data[[a]],
                Outlier       = iqr_flags(.data[[a]])) %>%
      ggplot(aes(x = MatrixType, y = Concentration)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(aes(color = Outlier, shape = Outlier),
                  width = 0.15, alpha = 0.85, size = 2.5) +
      scale_color_manual(values = c("FALSE" = "#555555", "TRUE" = "#e31a1c"),
                         na.value = "grey70") +
      scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 17), na.value = 4) +
      facet_wrap(~ Location) +
      labs(title    = paste("IQR Outliers —", a),
           x        = NULL,
           y        = "Concentration",
           color    = "IQR outlier",
           shape    = "IQR outlier") +
      theme_bw(base_size = 13)
  })

  output$mahal_tbl <- renderDT({
    md_results %>%
      head(10) %>%
      datatable(options = list(dom = "t", pageLength = 10), rownames = FALSE)
  })

  output$mahal_bar <- renderPlot({
    df_deg  <- length(analytes)
    co      <- c("p=0.05"  = qchisq(0.95,  df_deg),
                 "p=0.01"  = qchisq(0.99,  df_deg),
                 "p=0.001" = qchisq(0.999, df_deg))
    md_results %>%
      ggplot(aes(x = reorder(SampleID, Mahalanobis), y = Mahalanobis,
                 fill = MV_outlier)) +
      geom_col() +
      scale_fill_manual(values = c("FALSE" = "#4292c6", "TRUE" = "#e31a1c"),
                        labels = c("FALSE" = "No", "TRUE" = "Yes (p<0.001)")) +
      geom_hline(yintercept = co["p=0.05"],  color = "orange", linetype = "dashed", linewidth = 0.8) +
      geom_hline(yintercept = co["p=0.01"],  color = "red",    linetype = "dashed", linewidth = 0.8) +
      geom_hline(yintercept = co["p=0.001"], color = "purple", linetype = "dashed", linewidth = 0.8) +
      coord_flip() +
      labs(title    = "Mahalanobis Distance by Sample",
           subtitle = "Dashed lines: Chi² thresholds — orange=0.05, red=0.01, purple=0.001",
           x = NULL, y = "Mahalanobis distance", fill = "MV outlier") +
      theme_bw(base_size = 11)
  })

  output$fingerprint_plot <- renderPlot({
    least_id <- md_results %>% slice_min(Mahalanobis, n = 1) %>% pull(SampleID)
    most_id  <- md_results %>% slice_max(Mahalanobis, n = 1) %>% pull(SampleID)
    chem_rl2 %>%
      select(SampleID, all_of(analytes)) %>%
      filter(SampleID %in% c(least_id, most_id)) %>%
      pivot_longer(all_of(analytes), names_to = "Analyte", values_to = "Conc") %>%
      mutate(Type = ifelse(SampleID == least_id,
                           paste("Least abnormal:", least_id),
                           paste("Most abnormal:", most_id))) %>%
      ggplot(aes(x = Analyte, y = Conc, group = SampleID, color = Type)) +
      geom_line(linewidth = 1.1) +
      geom_point(size = 2.5) +
      labs(title = "Chemical Fingerprints: Least vs Most Abnormal (Mahalanobis)",
           x = NULL, y = "Concentration (RL/2-filled)", color = NULL) +
      theme_bw(base_size = 12) +
      theme(axis.text.x  = element_text(angle = 45, hjust = 1),
            legend.position = "bottom")
  })

  output$na_analyte_tbl <- renderDT({
    tibble(
      Analyte = analytes,
      N_NA    = map_int(analytes, ~ sum(is.na(chem_wide[[.x]]))),
      N_total = nrow(chem_wide)
    ) %>%
      mutate(NA_rate = round(N_NA / N_total, 3)) %>%
      arrange(desc(NA_rate)) %>%
      datatable(options = list(dom = "t", pageLength = 9), rownames = FALSE)
  })

  output$na_sample_tbl <- renderDT({
    chem_wide %>%
      transmute(SampleID,
                N_NA    = rowSums(is.na(across(all_of(analytes)))),
                NA_rate = round(N_NA / length(analytes), 2)) %>%
      arrange(desc(N_NA)) %>%
      head(10) %>%
      datatable(options = list(dom = "t"), rownames = FALSE)
  })

  output$na_heatmap <- renderPlot({
    chem_wide %>%
      select(SampleID, all_of(analytes)) %>%
      pivot_longer(-SampleID, names_to = "Analyte", values_to = "Conc") %>%
      mutate(Missing = is.na(Conc)) %>%
      ggplot(aes(x = Analyte, y = SampleID, fill = Missing)) +
      geom_tile() +
      scale_fill_manual(values = c("FALSE" = "#f0f0f0", "TRUE" = "#e31a1c")) +
      labs(title = "Missingness Map  (red = non-detect / NA)", x = NULL, y = NULL) +
      theme_bw(base_size = 12) +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  })

  output$imp_compare_tbl <- renderDT({
    a <- input$imp_analyte
    tibble(
      SampleID = chem_wide$SampleID,
      Observed = chem_wide[[a]],
      RL2      = round(chem_rl2[[a]], 3),
      LM       = round(chem_lm[[a]], 3),
      kNN      = round(chem_knn[[a]], 3),
      MICE     = round(chem_mice[[a]], 3),
      WasNA    = is.na(chem_wide[[a]])
    ) %>%
      filter(WasNA) %>%
      select(-WasNA) %>%
      datatable(options = list(dom = "t", pageLength = 10), rownames = FALSE)
  })

  output$imp_avp_plot <- renderPlot({
    method <- input$imp_method
    df     <- if (method == "lm") avp_lm else avp_knn
    label  <- if (method == "lm")
      "Least-squares: Actual vs Predicted (observed rows)"
    else
      "kNN: Actual vs Predicted (10% holdout validation)"

    df %>%
      ggplot(aes(x = Actual, y = Predicted)) +
      geom_point(alpha = 0.7) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
      facet_wrap(~ Analyte, scales = "free") +
      labs(title    = label,
           subtitle = "Dashed line = perfect prediction (y = x)",
           x = "Actual", y = "Predicted") +
      theme_bw(base_size = 11)
  })

  output$imp_resid_plot <- renderPlot({
    method <- input$imp_method
    df     <- if (method == "lm") avp_lm else avp_knn
    label  <- if (method == "lm") "Least-squares" else "kNN"

    df %>%
      ggplot(aes(x = Predicted, y = Residual)) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_point(alpha = 0.7) +
      facet_wrap(~ Analyte, scales = "free_x") +
      labs(title    = paste(label, "— Residual Diagnostics"),
           subtitle = "Random scatter around 0 = ideal; patterns indicate bias",
           x = "Predicted", y = "Residual (Actual − Predicted)") +
      theme_bw(base_size = 11)
  })

  output$dist_raw <- renderPlot({
    dist_plot(chem_rl2[[input$scale_analyte]], "Raw", "#4292c6")
  })
  output$dist_log <- renderPlot({
    dist_plot(log10(chem_rl2[[input$scale_analyte]] + 1e-6), "log₁₀", "#74c476")
  })
  output$dist_minmax <- renderPlot({
    x0 <- chem_rl2[[input$scale_analyte]]
    dist_plot((x0 - min(x0)) / (max(x0) - min(x0)), "Min-Max", "#fd8d3c")
  })
  output$dist_percent <- renderPlot({
    x0 <- chem_rl2[[input$scale_analyte]]
    dist_plot((x0 / sum(x0)) * 100, "Percent", "#9e9ac8")
  })

  output$shapiro_tbl <- renderDT({
    x0   <- chem_rl2[[input$scale_analyte]]
    xlog <- log10(x0 + 1e-6)
    xmm  <- (x0 - min(x0)) / (max(x0) - min(x0))
    xpct <- (x0 / sum(x0)) * 100

    tibble(
      Transform = c("Raw", "log10", "Min-Max", "Percent"),
      W         = map_dbl(list(x0, xlog, xmm, xpct),
                          ~ round(shapiro.test(.x)$statistic, 4)),
      p_value   = map_dbl(list(x0, xlog, xmm, xpct),
                          ~ round(shapiro.test(.x)$p.value, 4))
    ) %>%
      arrange(desc(p_value)) %>%
      datatable(options = list(dom = "t"), rownames = FALSE)
  })

  # ══════════════════════════════════════════════════════════
  # Chapter 5: PCA vs UMAP
  # ══════════════════════════════════════════════════════════

  output$ch5_content <- renderUI({
    if (is.null(nhanes)) {
      return(
        div(class = "alert alert-warning mt-4",
          h4("⚠️  NHANES.csv not found"),
          p("Place ", tags$code("NHANES.csv"), " in the same folder as ",
            tags$code("app.R"), " and restart the app."),
          tags$ul(
            tags$li("Analyte columns: 7 through 41"),
            tags$li("Must have an ", tags$code("Age"), " column for the color plot"),
            tags$li("No missing values in columns 7–41")
          )
        )
      )
    }
    tagList(
      plotOutput("ch5_side_by_side", height = "400px"),
      br(),
      plotOutput("ch5_umap_age", height = "420px")
    )
  })

  umap_result <- eventReactive(input$run_umap, {
    req(!is.null(nhanes))
    X <- nhanes[, 7:41]

    pca_res <- prcomp(X, center = FALSE, scale. = FALSE)
    pve     <- (pca_res$sdev^2) / sum(pca_res$sdev^2)

    cfg              <- umap.defaults
    cfg$n_neighbors  <- input$umap_k
    cfg$min_dist     <- input$umap_dist
    cfg$n_components <- 2
    umap_out <- umap(X, config = cfg)

    list(pca = pca_res, pve = pve, umap = umap_out)
  })

  output$ch5_side_by_side <- renderPlot({
    res  <- umap_result()
    pve1 <- round(100 * res$pve[1], 1)
    pve2 <- round(100 * res$pve[2], 1)

    pca_df  <- data.frame(D1 = res$pca$x[, 1],     D2 = res$pca$x[, 2])
    umap_df <- data.frame(D1 = res$umap$layout[, 1], D2 = res$umap$layout[, 2])

    p1 <- ggplot(pca_df, aes(D1, D2)) +
      geom_point(alpha = 0.8, size = 2) +
      theme_bw(base_size = 13) +
      labs(title = "PCA",
           x     = paste0("PC1 (", pve1, "%)"),
           y     = paste0("PC2 (", pve2, "%)"))

    p2 <- ggplot(umap_df, aes(D1, D2)) +
      geom_point(alpha = 0.8, size = 2) +
      theme_bw(base_size = 13) +
      labs(title    = "UMAP",
           subtitle = paste0("n_neighbors = ", input$umap_k,
                             "   min_dist = ", input$umap_dist),
           x = "UMAP 1", y = "UMAP 2")

    cowplot::plot_grid(p1, p2, ncol = 2)
  })

  output$ch5_umap_age <- renderPlot({
    res <- umap_result()
    df  <- data.frame(
      UMAP1 = res$umap$layout[, 1],
      UMAP2 = res$umap$layout[, 2],
      Age   = nhanes$Age
    )
    ggplot(df, aes(UMAP1, UMAP2, color = Age)) +
      geom_point(size = 2.8, alpha = 0.9) +
      scale_color_brewer(palette = "Set2") +
      theme_bw(base_size = 14) +
      labs(title    = "UMAP of NHANES Biomarkers — Colored by Age",
           subtitle = paste0("n_neighbors = ", input$umap_k,
                             "   min_dist = ", input$umap_dist),
           x = "UMAP 1", y = "UMAP 2", color = "Age Group")
  })
}

shinyApp(ui, server)
