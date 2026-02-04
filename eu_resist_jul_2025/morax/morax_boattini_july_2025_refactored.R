# morax_boattini_july_2025_refactored.R

library(readxl)
library(data.table)
library(tidyverse)
library(pheatmap)
options(scipen = 999)

# --- Config ---
CFG <- list(
  raw_xlsx = "../data/2.MORAX-EU_dataset.xlsx",
  lookup_csv = "../data/Look_up_MORAX_Complete_dataset.csv",
  output_dir = "../out",
  mic_range = c("Benzylpenicillin MIC","Doxycycline MIC"),
  zone_range = c("Benzylpenicillin inhib zone diam","Doxycycline inhib zone diam"),
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
  # ensure numeric thresholds
  num_cols <- intersect(c("EUCAST_mic_big","CLSI_mic_bigeq_m100","CASFM_mic_big","EUCAST_disc_lower","CLSI_disc_lower_m100","CASFM_disc_lower"), names(l))
  if (length(num_cols)>0) l <- l %>% mutate_at(vars(one_of(num_cols)), as.numeric)
  l
}

safe_read_raw <- function(path = CFG$raw_xlsx, sheet = "1. Complete dataset") {
  read_xlsx(path = path, sheet = sheet, col_types = "text", trim_ws = TRUE)
}

# sanitize MIC tokens (small map based on repetitive tokens seen in the original script)
mic_token_map <- c("≤1" = "0.1")

sanitize_mic <- function(x) {
  x <- as.character(x)
  x[x %in% names(mic_token_map)] <- mic_token_map[x[x %in% names(mic_token_map)]]
  as.numeric(x)
}

# --- MIC pipeline ---
process_mic <- function(raw, lookup) {
  mic_cols <- raw %>% select(contains(" MIC")) %>% names()
  if (length(mic_cols)==0) stop("No MIC columns found")

  mic_df <- raw %>% select(`Code event`, `Species identification`) %>%
    bind_cols(raw %>% select(all_of(mic_cols))) %>%
    gather(Abx, MIC, all_of(mic_cols)) %>%
    mutate(Abx = str_replace(Abx, " MIC", "")) %>%
    mutate(MIC_clean = sanitize_mic(MIC)) %>%
    left_join(lookup %>% select(-c(EUCAST_disc_lower, CLSI_disc_lower_m100, CASFM_disc_lower)), by = c("Species identification","Abx")) %>%
    mutate(EUCAST_mic_Resist = ifelse(!is.na(MIC_clean) & !is.na(EUCAST_mic_big) & MIC_clean > EUCAST_mic_big, 1, 0),
           CLIST_mic_m100_Resist = ifelse(!is.na(MIC_clean) & !is.na(CLSI_mic_bigeq_m100) & MIC_clean >= CLSI_mic_bigeq_m100, 1, 0),
           CASFM_mic_Resist = ifelse(!is.na(MIC_clean) & !is.na(CASFM_mic_big) & MIC_clean > CASFM_mic_big, 1, 0))
  mic_df
}

summarise_and_write <- function(df, resist_col, prefix) {
  summ <- df %>% filter(!is.na(!!sym(resist_col))) %>%
    group_by(`Species identification`, Abx) %>%
    summarise(n = n(), resistance_rate = mean(!!sym(resist_col))*100, .groups = "drop")
  fwrite(summ, file.path(CFG$output_dir, paste0(prefix, "_resistance_summary.csv")))
  # heatmap for top species by sample count
  top_species <- df %>% select(`Species identification`, `Code event`) %>% distinct() %>%
    group_by(`Species identification`) %>% count() %>% arrange(-n) %>% head(CFG$top_n) %>% pull(`Species identification`)
  heat_df <- summ %>% filter(`Species identification` %in% top_species) %>%
    pivot_wider(names_from = Abx, values_from = resistance_rate, values_fill = NA)
  mat <- as.matrix(heat_df[,-1])
  rownames(mat) <- heat_df[[1]]
  mat[is.na(mat)] <- -99
  display_mat <- round(mat,0)
  p <- pheatmap(mat, color = colorRampPalette(c("white","lightcyan","royalblue4"))(50),
                cluster_rows = TRUE, cluster_cols = TRUE, na_col = "white", display_numbers = display_mat,
                main = paste0(prefix, " Resistance [Top ",CFG$top_n,"]"))
  outfile <- file.path(CFG$output_dir, paste0("heatmap_", prefix, ".svg"))
  ggsave(outfile, plot = p, width = 8, height = 10)
  list(summary = summ, heatmap = outfile)
}

# --- ZONE pipeline ---
process_zone <- function(raw, lookup) {
  zone_cols <- raw %>% select(contains("inhib zone")) %>% names()
  if (length(zone_cols)==0) stop("No zone columns found")
  zone_df <- raw %>% select(`Code event`, `Species identification`) %>% bind_cols(raw %>% select(all_of(zone_cols))) %>%
    gather(Abx, ZONE, all_of(zone_cols)) %>%
    mutate(Abx = str_replace(Abx, " inhib zone diam", "")) %>%
    mutate(ZONE = as.numeric(ZONE)) %>%
    left_join(lookup %>% select(-c(EUCAST_mic_big, CLSI_mic_bigeq_m100, CASFM_mic_big)), by = c("Species identification","Abx")) %>%
    mutate(EUCAST_Diam_Resist = ifelse(!is.na(ZONE) & !is.na(EUCAST_disc_lower) & ZONE < EUCAST_disc_lower, 1, 0),
           CLIST_Diam_m100_Resist = ifelse(!is.na(ZONE) & !is.na(CLSI_disc_lower_m100) & ZONE < CLSI_disc_lower_m100, 1, 0),
           CASFM_Diam_Resist = ifelse(!is.na(ZONE) & !is.na(CASFM_disc_lower) & ZONE < CASFM_disc_lower, 1, 0))
  zone_df
}

# --- Main ---
main <- function() {
  message("Loading lookup and raw data...")
  lookup <- load_lookup(CFG$lookup_csv)
  raw <- safe_read_raw(CFG$raw_xlsx)

  message("Processing MICs")
  mic_df <- process_mic(raw, lookup)
  fwrite(mic_df, file.path(CFG$output_dir, "MIC_data_clean_Jul_16.csv"))
  mic_res <- summarise_and_write(mic_df, "EUCAST_mic_Resist", "mic_eucast")

  message("Processing Zone diameters")
  zone_df <- process_zone(raw, lookup)
  fwrite(zone_df, file.path(CFG$output_dir, "ZONE_data_clean_Jul_16.csv"))
  zone_res <- summarise_and_write(zone_df, "EUCAST_Diam_Resist", "zone_eucast")

  message("Writing resist counts and summaries")
  # write counts similar to original
  fwrite(mic_df %>% filter(!is.na(EUCAST_mic_Resist)) %>% group_by(`Species identification`, Abx, EUCAST_mic_Resist) %>% count() %>%
           pivot_wider(names_from = EUCAST_mic_Resist, values_from = n, values_fill = 0) %>% mutate(total = `0`+`1`) %>% mutate(perc = `1`/total),
         file.path(CFG$output_dir, "EUCAST_resist_counts_Jul_16.csv"))

  fwrite(zone_df %>% filter(!is.na(EUCAST_Diam_Resist)) %>% group_by(`Species identification`, Abx, EUCAST_Diam_Resist) %>% count() %>%
           pivot_wider(names_from = EUCAST_Diam_Resist, values_from = n, values_fill = 0) %>% mutate(total = `0`+`1`) %>% mutate(perc = `1`/total),
         file.path(CFG$output_dir, "EUCAST_resist_counts_zone_Jul_16.csv"))

  message("Refactor complete. Outputs in: ", CFG$output_dir)
}

if (interactive() || (!interactive() && identical(Sys.getenv("R_SCRIPT_RUNNING"), "1"))) main()
