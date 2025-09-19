# ============================ Dependencies ============================
# install.packages(c("tidyverse","janitor","ggrepel","scales"))
library(tidyverse)
library(janitor)
library(ggrepel)
library(scales)

# ============================ Configuration ============================
paths <- list(
  populations = "20131219.populations.tsv",    # 1KG populations.tsv
  metadata    = "pop_super_metadata.csv",      # no header: sample,pop,super_pop,gender
  bucket_map  = "bucket_map.csv",              # optional: two columns pop,bucket
  # trait_file  = "population_height.csv",       # or "population_bmi.csv"
  # score_file  = "scores/PGS002802_scores.txt"  # choose one PGS score file
  trait_file  = "population_bmi.csv",       # or "population_bmi.csv"
  score_file  = "scores/PGS000027_scores.txt"  # choose one PGS score file
)

# Output prefixes (the script will append _male / _female)
out_prefix_png <- "scatter_pgs_vs_bmi"
out_prefix_csv <- "points_pgs_vs_bmi"

# ============================ Helpers ============================
# Detect male/female columns automatically (handles Male_Avg / Female_Avg / Male / Female / AvgMale, etc.)
.pick_mf_cols <- function(nms){
  nn <- janitor::make_clean_names(nms)
  male_idx   <- which(stringr::str_detect(nn, "^male(_?avg|_?mean)?$|^avg_?male$|male"))
  female_idx <- which(stringr::str_detect(nn, "^female(_?avg|_?mean)?$|^avg_?female$|female"))
  list(
    male   = if (length(male_idx))   nms[male_idx[1]]   else NA_character_,
    female = if (length(female_idx)) nms[female_idx[1]] else NA_character_
  )
}

# Read a country-level trait table (two columns for male/female) → long format (male/female/both)
read_trait_long <- function(file, value_name = "trait_value",
                            code_col="Code", country_col="Country/Region"){
  stopifnot(file.exists(file))
  raw <- readr::read_csv(file, show_col_types = FALSE)
  mf <- .pick_mf_cols(names(raw))
  if (is.na(mf$male) || is.na(mf$female)) {
    stop("Could not find obvious male/female columns in: ", file, " (e.g., Male_Avg / Female_Avg).")
  }
  tb <- raw %>%
    clean_names() %>%
    rename(
      code       = all_of(janitor::make_clean_names(code_col)),
      country    = all_of(janitor::make_clean_names(country_col)),
      male_val   = all_of(janitor::make_clean_names(mf$male)),
      female_val = all_of(janitor::make_clean_names(mf$female))
    ) %>%
    mutate(across(c(male_val, female_val), ~ suppressWarnings(as.numeric(.))))
  
  long_sex <- tb %>%
    select(code, country, male_val, female_val) %>%
    pivot_longer(c(male_val, female_val), names_to = "sex", values_to = value_name) %>%
    mutate(sex = recode(sex, male_val="male", female_val="female"))
  
  both_mean <- tb %>%
    transmute(code, country, sex = "both",
              !!value_name := rowMeans(across(c(male_val, female_val)), na.rm = TRUE))
  
  bind_rows(long_sex, both_mean)
}

# Read a PGS score file (use the first non-`sample` column as score)
read_one_scores <- function(file) {
  stopifnot(file.exists(file))
  sr <- readr::read_csv(file, show_col_types = FALSE)
  score_cols <- setdiff(names(sr), "sample")
  stopifnot(length(score_cols) >= 1)
  nm <- score_cols[1]
  sr %>%
    select(sample, !!nm) %>%
    rename(score = !!nm) %>%
    mutate(score = suppressWarnings(as.numeric(score)))
}

# ============================ Load 1KG base data ============================
official_pop <- readr::read_tsv(paths$populations, show_col_types = FALSE) %>%
  clean_names() %>%
  rename(
    pop            = population_code,
    pop_name       = population_description,
    super_pop_name = super_population
  ) %>%
  mutate(
    super_pop_official = recode(super_pop_name,
                                "African"="AFR","Ad Mixed American"="AMR","East Asian"="EAS",
                                "European"="EUR","South Asian"="SAS", .default = NA_character_)
  ) %>%
  select(pop, pop_name, super_pop_name, super_pop_official)

metadata <- readr::read_csv(paths$metadata,
                            col_names = c("sample","pop","super_pop","gender"),
                            show_col_types = FALSE) %>%
  distinct() %>%
  mutate(gender_std = str_to_lower(str_trim(gender)),
         gender_std = recode(gender_std, "m"="male", "f"="female", .default = gender_std))

# pop → country mapping (keeps your custom mapping)
pop2country <- c(
  ACB="Barbados", ASW="United States of America", ESN="Nigeria",
  GWD="Gambia", LWK="Kenya", MSL="Sierra Leone", YRI="Nigeria",
  CLM="Colombia", MXL="Mexico", PEL="Peru", PUR="Puerto Rico",
  CDX="China", CHB="China", CHS="China", JPT="Japan", KHV="Viet Nam",
  CEU="United States of America",
  FIN="Finland", GBR="United Kingdom",
  IBS="Spain", TSI="Italy",
  BEB="Bangladesh", GIH="India", ITU="India",
  PJL="Pakistan", STU="Sri Lanka"
)

normalize_country <- function(x){
  dplyr::recode(x,
                "Vietnam"="Viet Nam","The Gambia"="Gambia",
                "USA"="United States of America","United States"="United States of America",
                "UK"="United Kingdom", .default = x)
}
missing_pops <- setdiff(unique(metadata$pop), names(pop2country))
if (length(missing_pops) > 0) {
  guess_tbl <- official_pop %>%
    filter(pop %in% missing_pops) %>%
    mutate(country_guess = str_trim(str_extract(pop_name, "[^,]+$")),
           country_guess = normalize_country(country_guess)) %>%
    select(pop, country_guess) %>% drop_na()
  if (nrow(guess_tbl) > 0) {
    add_map <- setNames(guess_tbl$country_guess, guess_tbl$pop)
    pop2country <- c(pop2country, add_map)
  }
}

# ============================ Read trait & scores ============================
trait_long <- read_trait_long(paths$trait_file, value_name = "trait_value")
scores     <- read_one_scores(paths$score_file)

# ============================ Palette ============================
pal <- c(AFR="#B39D00", AMR="#F8766D", EAS="#00BFC4", EUR="#00A9FF", SAS="#F564E3")

# ============================ Plot one sex helper ============================
plot_one_sex <- function(sex_label = c("male","female")) {
  sex_label <- match.arg(sex_label)
  
  # 1) Country-level trait for this sex
  trait_sel <- trait_long %>%
    filter(sex == sex_label) %>%
    select(country, trait_value)
  
  # 2) Filter samples by sex
  md <- metadata %>% filter(gender_std == sex_label)
  
  # 3) Merge samples, map country, optional bucket
  df <- md %>%
    left_join(official_pop, by="pop") %>%
    inner_join(scores, by="sample") %>%
    mutate(country = recode(pop, !!!as.list(pop2country))) %>%
    filter(!is.na(pop), !is.na(super_pop), !is.na(country))
  
  if (file.exists(paths$bucket_map)) {
    bucket_map <- readr::read_csv(paths$bucket_map, show_col_types = FALSE) %>%
      clean_names() %>% select(pop, bucket) %>% mutate(across(everything(), as.character))
    df <- df %>% left_join(bucket_map, by="pop") %>%
      mutate(bucket = if_else(is.na(bucket), pop, bucket))
  } else {
    df <- df %>% mutate(bucket = pop)
  }
  
  # 4) Mean PGS per bucket
  pgs_by_bucket <- df %>%
    group_by(pop, super_pop, bucket) %>%
    summarise(n=n(),
              mean_pgs = mean(score, na.rm = TRUE),
              sd_pgs   = sd(score, na.rm = TRUE),
              .groups="drop")
  
  # 5) Attach country-level trait to buckets
  bucket_info <- df %>% distinct(bucket, pop, super_pop, country)
  plot_df <- bucket_info %>%
    left_join(trait_sel, by="country") %>%
    left_join(pgs_by_bucket %>% select(bucket, super_pop, mean_pgs),
              by=c("bucket","super_pop")) %>%
    filter(!is.na(mean_pgs), !is.na(trait_value))
  
  if (nrow(plot_df) == 0) {
    warning("No points for sex = ", sex_label, "; skipped.")
    return(invisible(NULL))
  }
  
  # 6) Correlation & plot
  ct <- suppressWarnings(cor.test(plot_df$mean_pgs, plot_df$trait_value, method = "pearson"))
  r_txt <- sprintf("italic(r) == %.2f", unname(ct$estimate))
  p_txt <- sprintf("italic(p) == %s",
                   if (ct$p.value < 1e-4) format(ct$p.value, digits=1, scientific=TRUE)
                   else sprintf("%.3f", ct$p.value))
  
  p <- ggplot(plot_df, aes(x = mean_pgs, y = trait_value,
                           color = super_pop, label = bucket)) +
    geom_point(size = 3) +
    ggrepel::geom_text_repel(show.legend = FALSE, min.segment_length = 0) +
    scale_color_manual(values = pal, drop = FALSE) +
    geom_smooth(method = "lm", se = FALSE, color = "gray40", linetype = "dashed") +
    annotate("text", x = -Inf, y = Inf, label = r_txt, parse = TRUE,
             hjust = -0.1, vjust = 1.6, size = 5) +
    annotate("text", x = -Inf, y = Inf, label = p_txt, parse = TRUE,
             hjust = -0.1, vjust = 3.0, size = 5) +
    labs(
      title = paste0("PGS vs Trait (", str_to_title(sex_label), ")"),
      x = "PGS (mean by population)",
      y = "Trait value (country-level)",
      color = "super_pop"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      panel.background = element_rect(fill = "white", color = "black"),
      plot.background  = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "grey85"),
      panel.grid.minor = element_line(color = "grey92"),
      axis.line        = element_line(color = "black")
    )
  
  # 7) Save
  out_png <- sprintf("%s_%s.png", out_prefix_png, sex_label)
  out_csv <- sprintf("%s_%s.csv", out_prefix_csv, sex_label)
  ggsave(out_png, p, width = 6.5, height = 6, dpi = 300)
  readr::write_csv(plot_df, out_csv)
  message("Saved plot: ", out_png)
  message("Saved table: ", out_csv)
}

# ============================ Run: male + female ============================
plot_one_sex("male")
plot_one_sex("female")
