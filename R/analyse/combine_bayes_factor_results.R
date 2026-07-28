combine_bayes_factor_results <- function(
    wave1 = "outputs/tables/bayes_factor_results_wave1.csv",
    wave2 = "outputs/tables/bayes_factor_results_wave2.csv",
    wave3 = "outputs/tables/bayes_factor_results_wave3.csv",
    out = "outputs/tables/bayes_factor_results.csv") {
  
  w1 <- readr::read_csv(wave1,show_col_types = FALSE)
  w2 <- if (file.exists(wave2)) {readr::read_csv(wave2,show_col_types = FALSE)
  } else {w1[0, ]}
  w3 <- if (file.exists(wave3)) {readr::read_csv(wave3,show_col_types = FALSE)
  } else {w1[0, ]}
  
  readr::write_csv(dplyr::bind_rows(w1,w2,w3),out)
}