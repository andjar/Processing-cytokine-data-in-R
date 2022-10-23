import_data <- function(filename, get_all = FALSE) {
  if (is.null(filename)) {
    df <- data.table()
  } else {
    df <- fread(filename)

    if (!get_all) {
      df <- df[plate > 0 & batch > 0]
    }

    df[, unique_sample_ID := seq_len(nrow(df))]
    df[, bp := paste0("b", batch, "-p", plate)]
    df[, batch := factor(batch)]
    df[, plate := factor(plate)]
  }

  return(df)
}
