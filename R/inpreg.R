inpreg <- function(use_log = TRUE,
                   use_FI = TRUE,
                   filenames = list(
                       FI = "C:/Users/anderhja/Lokal_folder/git/fluorescens_alle.csv",
                       biorad = "C:/Users/anderhja/Lokal_folder/git/raw_alle.csv",
                       expected = "C:/Users/anderhja/Lokal_folder/git/fluorescens_expected.csv"
                     )
                   ) {
  obj <- InpregObject$new(use_log = use_log, use_FI = use_FI, filenames = filenames)
}
