# reads a single excel sheet obtained from Varioskan microplate reader
read_sheet <- function(data_file, sheet) {
  suppressMessages(
    {
      first_col <- read_excel(path = data_file,
                          sheet = sheet,
                          col_names = FALSE,
                          range = cell_cols("A"))[[1]]
    }
    )

  wl_range <- str_split(first_col[6], ",")[[1]]

  exp_type <- ifelse(grepl("Em", wl_range[1]), "emission", "excitation")
  fixed_wl <- as.numeric(gsub(".*?([0-9]+).*", "\\1", wl_range[2]))

  first_line <- which(first_col == "Wavelength")
  num_rows <- min(which(is.na(first_col[first_line:length(first_col)])))

  num_cols <- ncol(read_excel(path = data_file,
                              sheet = sheet,
                              col_names = TRUE,
                              range = cell_rows(first_line)))

  spectra <- read_excel(path = data_file,
                        sheet = sheet,
                        col_names = TRUE,
                        range = anchored("A11",
                                         dim = c(num_rows - 1, num_cols)))
  return(list(spectra = spectra, exp_type = exp_type, fixed_wl = fixed_wl))
}

# reads spectral data from multipage excel file
# obtained from Varioskan microplate reader
# returns data frame in long format
tidy_data <- function(spectra_lst) {
  data <- spectra_lst$spectra %>%
    pivot_longer(cols = -c("Wavelength"),
                 names_to = "cell",
                 values_to = "intensity") %>%
    rename(var_wavelength = Wavelength) %>%
    mutate(fixed_wavelength = spectra_lst$fixed_wl,
           exp_type = spectra_lst$exp_type,
           cell = str_extract(cell, "(?<=\\().+?(?=\\))"))

  return(data)
}

# reads samples data from excel file
# adding samples data to spectral data frame
# returns data frame in long format
read_plate <- function(file_name) {
  data_file <- here("data", file_name)

  data_sheets <- excel_sheets(path = data_file) %>%
    .[matches("^data_", vars = .)]

  samples <- read_excel(path = data_file,
                        sheet = "samples",
                        col_names = TRUE)

  data <- map(data_sheets, read_sheet, data_file = data_file) %>%
    map(tidy_data) %>%
    bind_rows() %>%
    left_join(samples, by = "cell")

  return(data)
}

# reads a csv file with emission spectrum from Ocean Optics spectrometer
# returns spectral data frame with two columns:
# wavelength (nm) and intensity (a.u.)
# note: decimal delimiter is ','
read_ocean_fluoro <- function(filename) {
  con <- file(description = filename, open = "r", blocking = TRUE)
  n <- 0
  chunk <- ""

  while (chunk != ">>>>>Begin Spectral Data<<<<<") {
    n <- n + 1
    chunk <- readLines(con = con, n = 1)
  }
  close(con)

  data <- read_delim(file = filename,
                     delim = "\t",
                     locale = locale(decimal_mark = ","),
                     skip = n,
                     col_names = c("wavelength", "intensity"),
                     show_col_types = FALSE)
  data$file <- tools::file_path_sans_ext(basename(filename))
  return(data)
}


# reads the csv file with absorbance spectrum from Ocean Optics spectrometer
# returns spectral data frame with two columns: wavelength (nm) and absorbance
# note: decimal delimiter is ','
read_ocean_abs <- function(filename) {
  con <- file(description = filename, open = "r", blocking = TRUE)
  n <- 0
  chunk <- ""

  while (chunk != ">>>>>Begin Spectral Data<<<<<") {
    n <- n + 1
    chunk <- readLines(con = con, n = 1)
  }
  close(con)

  data <- read_delim(file = filename,
                     delim = "\t",
                     locale = locale(decimal_mark = ","),
                     skip = n,
                     col_names = c("wavelength", "absorbance"),
                     show_col_types = FALSE)
  data$file <- tools::file_path_sans_ext(basename(filename))
  return(data)
}

# reads a sp file with emission spectrum from Perkin Elmer spectrometer
# returns spectral data frame with two columns:
# wavelength (nm) and intensity (a.u.)
read_sp_fluoro <- function(filename) {
  con <- file(description = filename, open = "r", blocking = TRUE)
  n <- 0
  chunk <- ""

  while (chunk != "#DATA") {
    n <- n + 1
    chunk <- readLines(con = con, n = 1)
  }
  close(con)

  data <- read_delim(file = filename,
                     delim = "\t",
                     locale = locale(decimal_mark = "."),
                     skip = n,
                     col_names = c("wavelength", "intensity"),
                     show_col_types = FALSE)
  data$intensity <- as.numeric(trimws(data$intensity))
  data$file <- tools::file_path_sans_ext(basename(filename))
  return(data)
}



bl_correction <- function(spectrum, wl_range) {
  offset <- spectrum %>%
    filter(wavelength > wl_range[1] & wavelength < wl_range[2]) %>%
    .[[2]] %>%
    mean()
  spectrum[, 2] <- spectrum[, 2] - offset
  return(spectrum)
}


smooth_runmed <- function(spectrum, k = 3) {
  spectrum[, 2] <- runmed(x = spectrum[[2]], k = k)
  return(spectrum)
}


integral <- function(spectrum, wl_range) {
  intensity <- spectrum %>%
    filter(wavelength > wl_range[1] & wavelength < wl_range[2]) %>%
    .[["intensity"]] %>%
    sum()
  return(tibble(integral = intensity / (wl_range[2] - wl_range[1])))
}

# a function to encode factorial variables into the [-1, 1] range
encode_var <- function(x) {
  coded <- 2 * (x - mean(x)) / (max(x) - min(x))
}
