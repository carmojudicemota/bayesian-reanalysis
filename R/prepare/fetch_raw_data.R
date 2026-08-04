sha256_file <- function(path) {
  if (requireNamespace("digest", quietly = TRUE)) {
    return(digest::digest(path, algo = "sha256", file = TRUE))
  }
  if (requireNamespace("openssl", quietly = TRUE)) {
    con <- file(path, "rb")
    on.exit(close(con))
    return(as.character(openssl::sha256(con)))
  }
  NA_character_
}

fetch_osf <- function(guid, dest_dir) {
  node <- osfr::osf_retrieve_node(guid)
  files <- osfr::osf_ls_files(node, n_max = Inf)
  if (nrow(files) == 0) {
    return("no files listed (project may be private or view-only)")
  }
  osfr::osf_download(files, path = dest_dir, recurse = TRUE, conflicts = "overwrite", progress = FALSE)
  sprintf("downloaded %d item(s) from osf.io/%s", nrow(files), guid)
}

fetch_url <- function(url, dest_dir) {
  dest <- file.path(dest_dir, basename(url))
  utils::download.file(url, dest, mode = "wb", quiet = TRUE)
  sprintf("downloaded %s", basename(url))
}

fetch_one <- function(source_type, locator, dest_dir) {
  switch(
    source_type,
    osf = fetch_osf(locator, dest_dir),
    url = fetch_url(locator, dest_dir),
    psycharchives = sprintf("manual: retrieve from https://doi.org/%s into %s", locator, dest_dir),
    manual = sprintf("manual: %s", locator),
    stop(sprintf("unknown source_type: %s", source_type))
  )
}

fetch_raw_data <- function(manifest_path = "config/raw_data_manifest.csv",
                           dest_root = "data/raw",
                           only = NULL) {
  if (!requireNamespace("osfr", quietly = TRUE)) {
    stop("Package 'osfr' is required. Install it and record it with renv::snapshot().")
  }
  manifest <- readr::read_csv(manifest_path, show_col_types = FALSE)
  if (!is.null(only)) {
    manifest <- dplyr::filter(manifest, .data$study_id %in% only)
  }
  results <- purrr::map_dfr(seq_len(nrow(manifest)), function(i) {
    row <- manifest[i, ]
    dest_dir <- file.path(dest_root, row$study_id)
    dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
    status <- tryCatch(
      fetch_one(row$source_type, row$locator, dest_dir),
      error = function(e) sprintf("error: %s", conditionMessage(e))
    )
    message(sprintf("[%s] %s", row$study_id, status))
    tibble::tibble(study_id = row$study_id, source_type = row$source_type,
                   locator = row$locator, status = status)
  })
  results
}

write_raw_checksums <- function(dest_root = "data/raw",
                                out_path = "config/raw_data_checksums.csv") {
  rel <- list.files(dest_root, recursive = TRUE)
  if (length(rel) == 0) {
    message("No files under ", dest_root)
    return(invisible(NULL))
  }
  full <- file.path(dest_root, rel)
  table <- tibble::tibble(
    path = rel,
    bytes = file.info(full)$size,
    sha256 = vapply(full, sha256_file, character(1))
  )
  readr::write_csv(table, out_path)
  message(sprintf("Wrote %d checksums to %s", nrow(table), out_path))
  invisible(table)
}

verify_raw_checksums <- function(checksums_path = "config/raw_data_checksums.csv",
                                 dest_root = "data/raw") {
  if (!file.exists(checksums_path)) {
    stop("Checksum file not found: ", checksums_path)
  }
  expected <- readr::read_csv(checksums_path, show_col_types = FALSE)
  full <- file.path(dest_root, expected$path)
  expected$actual <- vapply(full, function(p) if (file.exists(p)) sha256_file(p) else NA_character_, character(1))
  expected$ok <- !is.na(expected$actual) & expected$actual == expected$sha256
  bad <- dplyr::filter(expected, !.data$ok)
  if (nrow(bad) == 0) {
    message(sprintf("All %d files verified.", nrow(expected)))
  } else {
    message(sprintf("%d of %d files failed verification:", nrow(bad), nrow(expected)))
    for (p in bad$path) message("  ", p)
  }
  expected
}
