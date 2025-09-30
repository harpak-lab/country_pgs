# ============================ Dependencies ============================================
# install.packages(c("tidyverse","janitor","ggrepel","scales","broom"))
library(tidyverse)
library(janitor)
library(ggrepel)
library(scales)
library(broom)  # Add broom package for model result organization

# ============================ Configuration ============================================
paths <- list(
  gbd         = "../data/IHME-GBD_2021_DATA-cff1c199-1.csv",
  populations = "../data/20131219.populations.tsv",
  metadata    = "../data/pop_super_metadata.csv"
)

# PGS & Mapping file
paths$scores_dir <- "../scores"
paths$map        <- "../data/cause_pgs_map.csv"

# Configuration
target_year <- NA_integer_  # Target year, NA means use all years
min_sample_size <- 5       # Minimum sample size threshold (increased for Alzheimer's)
min_countries <- 2          # Minimum number of countries for modeling
alpha_level <- 0.05         # Significance level

# Output directories
out_dir_fig <- "../results/fig_by_cause_logit"
dir.create(out_dir_fig, showWarnings = FALSE, recursive = TRUE)
out_csv_all <- "../results/pop_meanPGS_vs_metric_ALLCAUSES_logit.csv"
out_model_summary <- "../results/logistic_model_summary.csv"

# ============================ Read and Process Data ========================================
# Read GBD data
gbd <- readr::read_csv(paths$gbd, show_col_types = FALSE) %>% clean_names()

# Standardize column names
colmap <- list(
  measure = intersect(names(gbd), c("measure_name","measure")),
  age     = intersect(names(gbd), c("age_name","age")),
  sex     = intersect(names(gbd), c("sex_name","sex")),
  metric  = intersect(names(gbd), c("metric_name","metric")),
  loc     = intersect(names(gbd), c("location_name","location")),
  year    = intersect(names(gbd), c("year")),
  value   = intersect(names(gbd), c("val","value")),
  cause   = intersect(names(gbd), c("cause_name","cause"))
)

# Verify required columns exist
required_cols <- c("measure", "age", "sex", "metric", "loc", "year", "value", "cause")
for(col in required_cols) {
  if(length(colmap[[col]]) == 0) {
    stop(paste("Required column", col, "not found in GBD data"))
  }
}

# Read population data and metadata
official_pop <- readr::read_tsv(paths$populations, show_col_types = FALSE) %>%
  clean_names() %>%
  rename(
    pop          = population_code,
    pop_name     = population_description,
    super_pop_name = super_population
  ) %>%
  mutate(
    super_pop_official = recode(super_pop_name,
                                "African"          = "AFR",
                                "Ad Mixed American"  = "AMR",
                                "East Asian"       = "EAS",
                                "European"         = "EUR",
                                "South Asian"      = "SAS",
                                .default = NA_character_)
  ) %>%
  select(pop, pop_name, super_pop_name, super_pop_official)

metadata <- readr::read_csv(paths$metadata,
                            col_names = c("sample","pop","super_pop","gender"),
                            show_col_types = FALSE) %>%
  distinct()

# ============================ Population to Country Mapping =================
pop2country <- c(
  ACB="Barbados", ASW="United States of America", ESN="Nigeria",
  GWD="Gambia", LWK="Kenya", MSL="Sierra Leone", YRI="Nigeria",
  CLM="Colombia", MXL="Mexico", PEL="Peru", PUR="Puerto Rico",
  CDX="China", CHB="China", CHS="China", JPT="Japan", KHV="Viet Nam",
  CEU="United States of America", FIN="Finland", GBR="United Kingdom",
  IBS="Spain", TSI="Italy",
  BEB="Bangladesh", GIH="India", ITU="India",
  PJL="Pakistan", STU="Sri Lanka"
)

normalize_country <- function(x){
  dplyr::recode(x,
                "Vietnam"="Viet Nam",
                "The Gambia"="Gambia",
                "USA"="United States of America",
                "United States"="United States of America",
                "UK"="United Kingdom",
                .default = x)
}

# Auto-guess missing population mappings
missing_pops <- setdiff(unique(metadata$pop), names(pop2country))
if (length(missing_pops) > 0) {
  guess_tbl <- official_pop %>%
    filter(pop %in% missing_pops) %>%
    mutate(country_guess = str_trim(str_extract(pop_name, "[^,]+$")),
           country_guess = normalize_country(country_guess)) %>%
    select(pop, country_guess) %>% 
    drop_na()
  
  if (nrow(guess_tbl) > 0) {
    add_map <- setNames(guess_tbl$country_guess, guess_tbl$pop)
    pop2country <- c(pop2country, add_map)
  }
}

# ============================ Helper Functions ========================================
# Helper function: automatically scan scores_dir and create mapping suggestions
create_mapping_suggestions <- function(cause_names, scores_dir) {
  if (!dir.exists(scores_dir)) {
    return(tibble(cause_name = cause_names, pgs_id = "", score_file = ""))
  }
  
  # Get all score files
  score_files <- list.files(scores_dir, pattern = "PGS\\d+_scores\\.txt$", full.names = FALSE)
  
  message("📂 Found ", length(score_files), " PGS score files in ", scores_dir)
  if (length(score_files) > 0) {
    message("   Example files: ", paste(head(score_files, 3), collapse = ", "))
  }
  
  suggestions <- tibble(cause_name = cause_names) %>%
    mutate(
      pgs_id = "",  # Users need to manually fill PGS ID
      score_file = "", # Or fill in filename
      available_pgs_ids = paste(extract_pgs_id_from_filename(score_files), collapse = ", ")
    )
  
  return(suggestions)
}

`%||%` <- function(x, y) if (is.null(x) || (length(x)==1 && is.na(x)) || (length(x)==1 && x=="")) y else x

guess_score_path <- function(pgs_id, scores_dir) {
  # Check input values
  if (is.null(pgs_id) || is.na(pgs_id) || pgs_id == "") {
    return("")
  }
  
  # Ensure scores_dir exists
  if (!dir.exists(scores_dir)) {
    warning("Scores directory does not exist: ", scores_dir)
    return("")
  }
  
  # Optimized search order for your file naming format
  candidates <- c(
    # Standard format: PGS[ID]_scores.txt (most common)
    file.path(scores_dir, paste0(pgs_id, "_scores.txt")),
    # Other possible formats
    file.path(scores_dir, paste0(pgs_id, "_scores.csv")),
    file.path(scores_dir, paste0(pgs_id, "scores.txt")),
    file.path(scores_dir, paste0(pgs_id, ".txt")),
    file.path(scores_dir, paste0(pgs_id, ".csv"))
  )
  
  existing <- candidates[file.exists(candidates)]
  if (length(existing) > 0) return(existing[1])
  return("")
}

# Helper function to extract PGS ID from filename
extract_pgs_id_from_filename <- function(filename) {
  # Extract "PGS002724" from filenames like "PGS002724_scores.txt"
  pgs_match <- str_extract(filename, "PGS\\d+")
  return(ifelse(is.na(pgs_match), "", pgs_match))
}

read_one_scores <- function(file) {
  if (!file.exists(file)) {
    warning("Score file does not exist: ", file)
    return(tibble())
  }
  
  sr <- readr::read_csv(file, show_col_types = FALSE)
  score_cols <- setdiff(names(sr), "sample")
  
  if (length(score_cols) == 0) {
    warning("No score columns found in file: ", file)
    return(tibble())
  }
  
  nm <- score_cols[1]
  sr %>%
    select(sample, !!nm) %>%
    rename(score = !!nm) %>%
    mutate(score = suppressWarnings(as.numeric(score))) %>%
    filter(!is.na(score))
}

base_filter <- function(d) {
  d %>%
    filter(
      str_to_lower(!!sym(colmap$age)) == "all ages",
      str_to_lower(!!sym(colmap$sex)) %in% c("both","both sexes"),
      str_detect(!!sym(colmap$metric), regex("rate", TRUE))
    ) %>%
    mutate(
      .value = suppressWarnings(as.numeric(!!sym(colmap$value))),
      .year  = suppressWarnings(as.integer(!!sym(colmap$year)))
    ) %>%
    filter(!is.na(.value), .value >= 0)
}

# ============================ Cause-PGS Mapping File ===========================
cause_map <- readr::read_csv(paths$map, show_col_types = FALSE) %>%
  clean_names()

# Ensure necessary columns exist
if (!"pgs_id" %in% names(cause_map)) {
  cause_map$pgs_id <- NA_character_
}
if (!"score_file" %in% names(cause_map)) {
  cause_map$score_file <- NA_character_
}

# Clean and validate mapping data
cause_map <- cause_map %>%
  mutate(
    cause_name = as.character(cause_name),
    pgs_id = case_when(
      is.na(pgs_id) ~ "",
      is.null(pgs_id) ~ "",
      TRUE ~ as.character(pgs_id)
    ),
    score_file = case_when(
      is.na(score_file) ~ "",
      is.null(score_file) ~ "",
      TRUE ~ as.character(score_file)
    )
  ) %>%
  filter(
    !is.na(cause_name), 
    nchar(trimws(cause_name)) > 0,
    # At least one of pgs_id or score_file must be provided
    (nchar(trimws(pgs_id)) > 0 | nchar(trimws(score_file)) > 0)
  ) %>%
  distinct()

if (nrow(cause_map) == 0) {
  stop("No valid mappings in mapping file. Please ensure each row has cause_name and (pgs_id or score_file).")
}

message("✅ Successfully loaded ", nrow(cause_map), " disease-PGS mappings")


# ============================ Logistic Regression =========================
run_logistic_analysis <- function(cname, measure_key = "prevalence") {
  cat("Analyzing disease:", cname, "- measure:", measure_key, "\n")
  
  # 1) Locate PGS score file
  row <- cause_map %>% 
    filter(str_to_lower(cause_name) == str_to_lower(cname)) %>% 
    slice(1)
  
  if (nrow(row) == 0) {
    message("⚠️ Disease not found in mapping file: ", cname)
    return(tibble())
  }
  
  sfile <- ""
  
  # Try to get from score_file column
  if ("score_file" %in% names(row)) {
    score_file_val <- row$score_file[1]
    if (!is.na(score_file_val) && !is.null(score_file_val) && score_file_val != "") {
      score_file_val <- as.character(score_file_val)
      
      # If it's an absolute path or relative path (contains path separator)
      if (str_detect(score_file_val, "[/\\\\]") && file.exists(score_file_val)) {
        sfile <- score_file_val
      } 
      # If it's just a filename, concatenate with scores_dir
      else {
        candidate_path <- file.path(paths$scores_dir, score_file_val)
        if (file.exists(candidate_path)) {
          sfile <- candidate_path
        }
      }
    }
  }
  
  # If file not found yet, try using pgs_id to guess
  if (sfile == "") {
    pgs_id_val <- ""
    if ("pgs_id" %in% names(row)) {
      pgs_id_temp <- row$pgs_id[1]
      if (!is.na(pgs_id_temp) && !is.null(pgs_id_temp)) {
        pgs_id_val <- as.character(pgs_id_temp)
      }
    }
    
    if (pgs_id_val != "") {
      guessed_path <- guess_score_path(pgs_id_val, paths$scores_dir)
      if (guessed_path != "" && file.exists(guessed_path)) {
        sfile <- guessed_path
      }
    }
  }
  
  # Final check if file exists
  if (sfile == "" || !file.exists(sfile)) {
    # List files in scores_dir for debugging
    if (dir.exists(paths$scores_dir)) {
      available_files <- list.files(paths$scores_dir, pattern = "\\.(txt|csv)$", full.names = FALSE)
      message("⚠️ PGS score file not found (disease: ", cname, ")")
      message("   Search path: ", paths$scores_dir)
      message("   Available files: ", paste(head(available_files, 5), collapse = ", "), 
              if(length(available_files) > 5) "..." else "")
    } else {
      message("⚠️ Score file directory does not exist: ", paths$scores_dir)
    }
    return(tibble())
  }
  
  # 2) Read scores and merge with metadata
  scores <- read_one_scores(sfile)
  if (nrow(scores) == 0) {
    message("⚠️ Score file is empty or invalid: ", sfile)
    return(tibble())
  }
  
  df_local <- metadata %>%
    left_join(official_pop, by = "pop") %>%
    inner_join(scores, by = "sample") %>%
    mutate(country = dplyr::recode(pop, !!!as.list(pop2country))) %>%
    filter(!is.na(pop), !is.na(super_pop), !is.na(country), !is.na(score)) %>%
    mutate(
      cause_name = cname,
      pgs_id = if (is.na(row$pgs_id) || row$pgs_id == "") basename(sfile) else row$pgs_id
    )
  
  if (nrow(df_local) == 0) {
    message("⚠️ No data after merging: ", cname)
    return(tibble())
  }
  
  # 3) Get GBD data
  gbd_sub <- gbd %>%
    filter(
      str_to_lower(!!sym(colmap$cause)) == str_to_lower(cname),
      str_detect(!!sym(colmap$measure), regex(measure_key, TRUE))
    ) %>%
    base_filter()
  
  if (!is.na(target_year)) {
    gbd_sub <- gbd_sub %>% filter(.year == target_year)
  }
  
  if (nrow(gbd_sub) == 0) {
    message("⚠️ No GBD data: ", cname, " - ", measure_key)
    return(tibble())
  }
  
  # Country-level prevalence data
  gbd_country_local <- gbd_sub %>%
    group_by(country = !!sym(colmap$loc)) %>%
    summarise(
      y_metric = mean(.value, na.rm = TRUE),
      n_years = n_distinct(.year),
      .groups = "drop"
    )
  
  # 4) Calculate population-level PGS statistics
  pgs_by_pop_local <- df_local %>%
    group_by(country, cause_name, pgs_id, pop, super_pop) %>%
    summarise(
      n_samples = n(),
      mean_pgs  = mean(score, na.rm = TRUE),
      sd_pgs    = sd(score, na.rm = TRUE),
      median_pgs = median(score, na.rm = TRUE),
      .groups   = "drop"
    ) %>%
    filter(n_samples >= min_sample_size)
  
  if (nrow(pgs_by_pop_local) == 0) {
    message("⚠️ Insufficient sample size (< ", min_sample_size, "): ", cname)
    return(tibble())
  }
  
  # 5) Country-level data aggregation
  country_tbl <- pgs_by_pop_local %>%
    group_by(country) %>%
    summarise(
      N = sum(n_samples),
      mean_pgs_country = weighted.mean(mean_pgs, w = n_samples, na.rm = TRUE),
      super_pop_country = names(sort(tapply(n_samples, super_pop, sum), decreasing = TRUE))[1],
      n_populations = n_distinct(pop),
      .groups = "drop"
    )
  
  # 6) Merge with prevalence data
  eps <- 1e-9
  combo_country <- country_tbl %>%
    left_join(gbd_country_local, by = "country") %>%
    filter(!is.na(y_metric), N > 0) %>%
    mutate(
      # Convert to probability (0,1)
      p = pmin(pmax(y_metric / 1e5, eps), 1 - eps),
      # Binomial distribution counts
      cases = pmin(pmax(floor(p * N), 0L), N),
      fails = N - cases,
      # Add quality control flag
      is_reliable = N >= min_sample_size & cases > 0 & fails > 0
    ) %>%
    filter(is_reliable)
  
  if (nrow(combo_country) < min_countries) {
    message("⚠️ Insufficient countries for modeling (< ", min_countries, "): ", cname, 
            " (found: ", nrow(combo_country), " countries)")
    
    # Special handling for diseases with limited data points (like Alzheimer's)
    if (nrow(combo_country) > 0) {
      message("   Available countries: ", paste(combo_country$country, collapse = ", "))
      message("   Sample sizes: ", paste(combo_country$N, collapse = ", "))
    }
    return(tibble())
  }
  
  # 7) PGS standardization (based on country means)
  mu_pgs <- mean(combo_country$mean_pgs_country, na.rm = TRUE)
  sd_pgs <- sd(combo_country$mean_pgs_country, na.rm = TRUE)
  
  if (!is.finite(sd_pgs) || sd_pgs <= 0) {
    message("⚠️ Invalid PGS standard deviation: ", cname)
    return(tibble())
  }
  
  combo_country <- combo_country %>%
    mutate(pgs_z_country = (mean_pgs_country - mu_pgs) / sd_pgs)
  
  pgs_by_pop_local <- pgs_by_pop_local %>%
    mutate(pgs_z_pop = (mean_pgs - mu_pgs) / sd_pgs)
  
  # 8) Logistic regression modeling
  tryCatch({
    fit <- glm(cbind(cases, fails) ~ pgs_z_country, 
               family = binomial(link = "logit"), 
               data = combo_country)
    
    # Model diagnostics
    model_summary <- broom::tidy(fit, conf.int = TRUE, conf.level = 1 - alpha_level)
    model_glance <- broom::glance(fit)
    
    # Extract PGS coefficient
    pgs_coef <- model_summary %>% filter(term == "pgs_z_country")
    
    if (nrow(pgs_coef) == 0) {
      message("⚠️ Cannot extract PGS coefficient: ", cname)
      return(tibble())
    }
    
    or <- exp(pgs_coef$estimate)
    or_ci_lo <- exp(pgs_coef$conf.low)
    or_ci_hi <- exp(pgs_coef$conf.high)
    pval <- pgs_coef$p.value
    
    # Model goodness of fit
    pseudo_r2 <- 1 - (model_glance$deviance / model_glance$null.deviance)
    aic <- model_glance$AIC
    
    message(sprintf("✅ %s | OR = %.3f [%.3f, %.3f]; p = %.3g; Pseudo-R² = %.3f; N countries = %d",
                    cname, or, or_ci_lo, or_ci_hi, pval, pseudo_r2, nrow(combo_country)))
    
    # 9) Generate prediction data for plotting
    pgs_range_pop <- range(pgs_by_pop_local$pgs_z_pop, na.rm = TRUE)
    pgs_range_country <- range(combo_country$pgs_z_country, na.rm = TRUE)
    pgs_range_all <- c(min(pgs_range_country[1], pgs_range_pop[1]),
                       max(pgs_range_country[2], pgs_range_pop[2]))
    newd <- tibble(
      pgs_z_country = seq(pgs_range_all[1], pgs_range_all[2], length.out = 200)
    )
    
    pred_results <- predict(fit, newdata = newd, type = "response", se.fit = TRUE)
    newd$pred_p <- pred_results$fit
    newd$pred_se <- pred_results$se.fit
    newd$pred_ci_lo <- plogis(qlogis(newd$pred_p) - 1.96 * newd$pred_se)
    newd$pred_ci_hi <- plogis(qlogis(newd$pred_p) + 1.96 * newd$pred_se)
    
    # 10) Create plotting data
    plot_df_pop <- pgs_by_pop_local %>%
      select(country, pop, super_pop, pgs_z_pop, n_samples) %>%
      left_join(
        gbd_country_local %>% select(country, y_metric),
        by = "country"
      ) %>%
      filter(!is.na(y_metric))%>%
      mutate(
        p_obs = y_metric / 1e5,
        point_type = "population"
      )
    
    plot_df_country <- combo_country %>%
      select(country, pgs_z_country, p, N, super_pop_country) %>%
      rename(
        pgs_z = pgs_z_country,
        p_obs = p,
        n_samples = N,
        super_pop = super_pop_country
      ) %>%
      mutate(
        pop = country,
        point_type = "country"
      )
    
    # 11) Create scatter plot with fitted curve
    pal <- c(AFR="#B39D00", AMR="#F8766D", EAS="#00BFC4", EUR="#00A9FF", SAS="#F564E3")
    
    p_plot <- ggplot() +

      geom_line(data = newd, 
                aes(x = pgs_z_country, y = pred_p),
                linewidth = 1.2, color = "black") +
      geom_point(data = plot_df_pop,
                 aes(x = pgs_z_pop, y = p_obs, color = super_pop, size = n_samples), 
                 alpha = 0.7, shape = 16) +
      # geom_point(data = plot_df_country,
      #            aes(x = pgs_z, y = p_obs, color = super_pop, size = n_samples),
      #            shape = 18, alpha = 0.9, stroke = 1.5) +
      geom_text_repel(data = plot_df_pop,
                      aes(x = pgs_z_pop, y = p_obs, label = pop, color = super_pop), 
                      show.legend = FALSE, size = 5, 
                      min.segment.length = 0, max.overlaps = 20) +
      
      scale_color_manual(values = pal, drop = FALSE, name = "Super Population") +
      scale_size_continuous(name = "Sample Size", range = c(2, 8)) +
      scale_y_continuous(labels = scales::percent_format(), 
                         limits = c(0, max(c(plot_df_pop$p_obs, plot_df_country$p_obs)) * 1.1)) +
      labs(
        title = paste0(cname, " | ", str_to_title(measure_key)),
        subtitle = sprintf(
          "Country-level Logistic Model: OR per 1-SD PGS = %.2f [%.2f, %.2f]; p = %s; Pseudo-R² = %.3f\nN = %d countries | Circles = populations, Diamonds = country aggregates",
          or, or_ci_lo, or_ci_hi,
          if (pval < 1e-4) format(pval, digits = 1, scientific = TRUE) else sprintf("%.4f", pval),
          pseudo_r2, nrow(combo_country)
        ),
        x = "Polygenic Score (Z-score)",
        y = paste(str_to_title(measure_key), "(Probability)"),
        caption = sprintf("Total populations = %d | Total samples = %d | Model fitted on country-level aggregated data",
                          nrow(plot_df_pop), sum(plot_df_pop$n_samples))
      ) +
      theme_minimal(base_size = 12) +
      theme(
        panel.background = element_rect(fill = "white", color = "black"),
        plot.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
        panel.grid.minor = element_line(color = "grey95", linewidth = 0.3),
        axis.line = element_line(color = "black"),
        legend.position = "bottom",
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 10, lineheight = 1.2)
      )
    
    # Save plot
    clean_name <- janitor::make_clean_names(cname)
    out_png <- file.path(out_dir_fig, 
                         paste0("logistic_", clean_name, "_", measure_key, ".png"))
    ggsave(out_png, p_plot, width = 11, height = 8, dpi = 300)
    
    # 12) Return results
    result <- tibble(
      cause_name = cname,
      measure = measure_key,
      n_countries = nrow(combo_country),
      n_populations = nrow(plot_df_pop),
      total_samples = sum(combo_country$N),
      or_per_sd_pgs = or,
      or_ci_lo = or_ci_lo,
      or_ci_hi = or_ci_hi,
      pval_pgs = pval,
      coefficient = pgs_coef$estimate,
      std_error = pgs_coef$std.error,
      pseudo_r2 = pseudo_r2,
      aic = aic,
      converged = fit$converged,
      pgs_mean = mu_pgs,
      pgs_sd = sd_pgs,
      model_type = "binomial_logistic_country_level",
      significant = pval < alpha_level,
      plot_file = basename(out_png),
      data_quality = case_when(
        nrow(combo_country) >= 8 ~ "Good",
        nrow(combo_country) >= min_countries ~ "Acceptable", 
        TRUE ~ "Limited"
      ),
      within_country_pgs_var = var(pgs_by_pop_local$pgs_z_pop, na.rm = TRUE),
      between_country_pgs_var = var(combo_country$pgs_z_country, na.rm = TRUE)
    )
    
    return(result)
    
  }, error = function(e) {
    message("⚠️ Model fitting failed: ", cname, " - ", e$message)
    return(tibble())
  })
}
# ============================ Batch Analysis =========================
message("Starting batch logistic regression analysis...")

all_results <- list()
skipped_diseases <- character()

for (cname in cause_map$cause_name) {
    print(cname)
    result <- run_logistic_analysis(cname, "prevalence")
    if (nrow(result) > 0) {
      all_results[[length(all_results) + 1]] <- result
    } else {
      skipped_diseases <- c(skipped_diseases, cname)
  }
}

# Combine all results
if (length(all_results) > 0) {
  final_results <- bind_rows(all_results)
  
  # Save detailed results
  readr::write_csv(final_results, out_csv_all)
  
  # Create simplified summary table
  summary_table <- final_results %>%
    arrange(pval_pgs) %>%
    select(cause_name, or_per_sd_pgs, or_ci_lo, or_ci_hi, pval_pgs, 
           pseudo_r2, n_countries, total_samples, significant, data_quality) %>%
    mutate(
      or_ci = sprintf("%.3f [%.3f, %.3f]", or_per_sd_pgs, or_ci_lo, or_ci_hi),
      p_formatted = case_when(
        pval_pgs < 0.001 ~ "<0.001",
        pval_pgs < 0.01 ~ sprintf("%.3f", pval_pgs),
        TRUE ~ sprintf("%.3f", pval_pgs)
      )
    )
  
  readr::write_csv(summary_table, out_model_summary)
  
  # Output summary statistics
  n_significant <- sum(final_results$significant, na.rm = TRUE)
  n_total <- nrow(final_results)
  n_good_quality <- sum(final_results$data_quality == "Good", na.rm = TRUE)
  n_acceptable_quality <- sum(final_results$data_quality == "Acceptable", na.rm = TRUE)
  
  message(sprintf("\n=== Analysis Complete ==="))
  message(sprintf("Total diseases analyzed: %d", n_total))
  message(sprintf("Significant associations (p < %.3f): %d (%.1f%%)", 
                  alpha_level, n_significant, n_significant/n_total*100))
  message(sprintf("Data quality - Good: %d, Acceptable: %d, Limited: %d", 
                  n_good_quality, n_acceptable_quality, 
                  n_total - n_good_quality - n_acceptable_quality))
  
  if (length(skipped_diseases) > 0) {
    message(sprintf("Diseases skipped due to insufficient data: %d", length(skipped_diseases)))
    message("   Skipped diseases: ", paste(head(skipped_diseases, 5), collapse = ", "))
    if (length(skipped_diseases) > 5) message("   (and ", length(skipped_diseases) - 5, " others)")
    
    # Special note about Alzheimer's and similar conditions
    alzheimer_related <- skipped_diseases[str_detect(str_to_lower(skipped_diseases), 
                                                   "alzheimer|dementia|cognitive")]
    if (length(alzheimer_related) > 0) {
      message("\n📝 Note: Alzheimer's and related cognitive disorders often have limited")
      message("   country-level prevalence data due to underdiagnosis and reporting")
      message("   differences across populations. Consider:")
      message("   - Using incidence data if available")
      message("   - Lowering min_sample_size threshold")
      message("   - Combining related phenotypes")
    }
  }
  
  message(sprintf("\nResults saved to: %s", out_csv_all))
  message(sprintf("Summary table saved to: %s", out_model_summary))
  message(sprintf("Plots saved to: %s", out_dir_fig))
  
  # Display top significant results
  if (n_significant > 0) {
    message("\n🔍 Top significant associations:")
    top_results <- summary_table %>%
      filter(significant) %>%
      arrange(pval_pgs) %>%
      head(5) %>%
      select(cause_name, or_ci, p_formatted, n_countries, data_quality)
    
    print(top_results, n = Inf)
  }
  
} else {
  message("⚠️ No diseases were successfully analyzed")
  message("Common issues:")
  message("  - PGS score files not found in specified directory")
  message("  - Insufficient sample sizes (< ", min_sample_size, " per population)")
  message("  - Insufficient countries (< ", min_countries, " for modeling)")
  message("  - Missing GBD prevalence data for specified diseases")
  message("\nTroubleshooting suggestions:")
  message("  1. Check that score files exist in: ", paths$scores_dir)
  message("  2. Verify mapping file contains correct PGS IDs or file names")
  message("  3. Consider lowering min_sample_size (currently ", min_sample_size, ")")
  message("  4. Consider lowering min_countries (currently ", min_countries, ")")
}

