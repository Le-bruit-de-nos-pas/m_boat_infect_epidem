# minofeu_boattini_july_2025_refactored.R

library(readxl)
library(data.table)
library(tidyverse)
library(pheatmap)
options(scipen = 999)

# --- Config ---
CFG <- list(
  raw_xlsx = "../data/2.MINOFEu_dataset_matteo_only.xlsx",
  lookup_csv = "../data/Look_up_MinoFeu_Complete_dataset.csv",
  output_dir = "../out",
  mic_pattern = " MIC",
  zone_pattern = " inhib zone",
  top_n = 50,
  min_samples = 10,
  seed = 123
)

if (!dir.exists(CFG$output_dir)) dir.create(CFG$output_dir, recursive = TRUE)

# --- Helpers ---
load_lookup <- function(path = CFG$lookup_csv) {
  if (!file.exists(path)) stop("Lookup CSV not found: ", path)
  l <- fread(path, colClasses = "character") %>%
    mutate(`Species identification` = str_replace_all(`Species identification`, "�", " ")) %>%
    mutate(`Species identification` = str_trim(`Species identification`)) %>% distinct()
  num_cols <- intersect(c("EUCAST_mic_big","CLSI_mic_bigeq_m100","CLSI_mic_bigeq_m45","CASFM_mic_big","EUCAST_disc_lower","CLSI_disc_lower_m100","CLSI_disc_lower_m45","CASFM_disc_lower"), names(l))
  if (length(num_cols)>0) l <- l %>% mutate_at(vars(one_of(num_cols)), as.numeric)
  l
}

safe_read_raw <- function(path = CFG$raw_xlsx, sheet = "3. Complete dataset") {
  read_xlsx(path = path, sheet = sheet, col_types = "text", trim_ws = TRUE)
}

# minimal sanitizer (original file appears numeric-clean for MICs but keep a hook)
sanitize_numeric <- function(x) {
  x <- as.character(x)
  x[x == ""] <- NA
  suppressWarnings(as.numeric(x))
}

# --- MIC processing ---
process_mic <- function(raw, lookup) {
  mic_cols <- raw %>% select(contains(CFG$mic_pattern)) %>% names()
  if (length(mic_cols) == 0) stop("No MIC columns detected")
  mic_df <- raw %>% select(`Code event`, `Species identification`) %>%
    bind_cols(raw %>% select(all_of(mic_cols))) %>%
    gather(Abx, MIC, all_of(mic_cols)) %>%
    mutate(Abx = str_replace(Abx, CFG$mic_pattern, "")) %>%
    mutate(MIC_clean = sanitize_numeric(MIC)) %>%
    left_join(lookup %>% select(-c(EUCAST_disc_lower, CLSI_disc_lower_m100, CLSI_disc_lower_m45, CASFM_disc_lower)), by = c("Species identification","Abx")) %>%
    mutate(EUCAST_mic_Resist = ifelse(!is.na(MIC_clean) & !is.na(EUCAST_mic_big) & MIC_clean > EUCAST_mic_big, 1, 0),
           CLIST_mic_m100_Resist = ifelse(!is.na(MIC_clean) & !is.na(CLSI_mic_bigeq_m100) & MIC_clean >= CLSI_mic_bigeq_m100, 1, 0),
           CLIST_mic_m45_Resist = ifelse(!is.na(MIC_clean) & !is.na(CLSI_mic_bigeq_m45) & MIC_clean >= CLSI_mic_bigeq_m45, 1, 0),
           CASFM_mic_Resist = ifelse(!is.na(MIC_clean) & !is.na(CASFM_mic_big) & MIC_clean > CASFM_mic_big, 1, 0))
  mic_df
}

# --- ZONE processing ---
process_zone <- function(raw, lookup) {
  zone_cols <- raw %>% select(contains(CFG$zone_pattern)) %>% names()
  if (length(zone_cols)==0) stop("No zone columns found")
  zone_df <- raw %>% select(`Code event`, `Species identification`) %>% bind_cols(raw %>% select(all_of(zone_cols))) %>%
    gather(Abx, ZONE, all_of(zone_cols)) %>%
    mutate(Abx = str_replace(Abx, " inhib zone diam", "")) %>%
    mutate(ZONE_clean = sanitize_numeric(ZONE)) %>%
    left_join(lookup %>% select(-c(EUCAST_mic_big, CLSI_mic_bigeq_m100, CLSI_mic_bigeq_m45, CASFM_mic_big)), by = c("Species identification","Abx")) %>%
    mutate(EUCAST_Diam_Resist = ifelse(!is.na(ZONE_clean) & !is.na(EUCAST_disc_lower) & ZONE_clean < EUCAST_disc_lower, 1, 0),
           CLIST_Diam_m100_Resist = ifelse(!is.na(ZONE_clean) & !is.na(CLSI_disc_lower_m100) & ZONE_clean < CLSI_disc_lower_m100, 1, 0),
           CLIST_Diam_m45_Resist = ifelse(!is.na(ZONE_clean) & !is.na(CLSI_disc_lower_m45) & ZONE_clean < CLSI_disc_lower_m45, 1, 0),
           CASFM_Diam_Resist = ifelse(!is.na(ZONE_clean) & !is.na(CASFM_disc_lower) & ZONE_clean < CASFM_disc_lower, 1, 0))
  zone_df
}

# --- Summaries & plots ---
summarize_resist <- function(df, resist_col, prefix) {
  summ <- df %>% filter(!is.na(!!sym(resist_col))) %>%
    group_by(`Species identification`, Abx) %>%
    summarise(n_samples = n(), resistance_rate = mean(!!sym(resist_col), na.rm = TRUE)*100, .groups = "drop")
  fwrite(summ, file.path(CFG$output_dir, paste0(prefix, "_resistance_summary.csv")))

  top_species <- df %>% select(`Species identification`, `Code event`) %>% distinct() %>%
    group_by(`Species identification`) %>% count() %>% arrange(-n) %>% head(CFG$top_n) %>% pull(`Species identification`)

  heat_df <- summ %>% filter(`Species identification` %in% top_species) %>%
    pivot_wider(names_from = Abx, values_from = resistance_rate, values_fill = NA)
  mat <- as.matrix(heat_df[,-1])
  rownames(mat) <- heat_df[[1]]
  mat[is.na(mat)] <- -99
  display_mat <- round(mat,0)
  p <- pheatmap(mat, color = colorRampPalette(c("white","lightcyan","royalblue4"))(50), cluster_rows = TRUE, cluster_cols = TRUE, na_col = "white", display_numbers = display_mat, main = paste0(prefix, " resistance [Top ", CFG$top_n, "]"))
  out <- file.path(CFG$output_dir, paste0("heatmap_", prefix, ".svg"))
  ggsave(out, plot = p, width = 7, height = 10)
  list(summary = summ, heatmap = out)
}

# --- Main ---
main <- function() {
  message("Loading lookup and raw dataset...")
  lookup <- load_lookup(CFG$lookup_csv)
  raw <- safe_read_raw(CFG$raw_xlsx)

  message("Processing MICs...")
  mic_df <- process_mic(raw, lookup)
  fwrite(mic_df, file.path(CFG$output_dir, "MIC_data_clean_minofeu.csv"))
  mic_res <- summarize_resist(mic_df, "EUCAST_mic_Resist", "mic_eucast")

  message("Processing zones...")
  zone_df <- process_zone(raw, lookup)
  fwrite(zone_df, file.path(CFG$output_dir, "ZONE_data_clean_minofeu.csv"))
  zone_res <- summarize_resist(zone_df, "EUCAST_Diam_Resist", "zone_eucast")

  message("Writing counts (EUCAST/CLSI/CASFM)...")
  fwrite(mic_df %>% filter(!is.na(EUCAST_mic_Resist)) %>% group_by(`Species identification`, Abx, EUCAST_mic_Resist) %>% count() %>% pivot_wider(names_from = EUCAST_mic_Resist, values_from = n, values_fill = 0) %>% mutate(total = `0`+`1`, perc = `1`/total), file.path(CFG$output_dir, "EUCAST_resist_counts_minofeu.csv"))

  fwrite(zone_df %>% filter(!is.na(EUCAST_Diam_Resist)) %>% group_by(`Species identification`, Abx, EUCAST_Diam_Resist) %>% count() %>% pivot_wider(names_from = EUCAST_Diam_Resist, values_from = n, values_fill = 0) %>% mutate(total = `0`+`1`, perc = `1`/total), file.path(CFG$output_dir, "EUCAST_resist_counts_zone_minofeu.csv"))

  message("Refactor complete. Outputs in:", CFG$output_dir)
}

if (interactive() || (!interactive() && identical(Sys.getenv("R_SCRIPT_RUNNING"), "1"))) main()
