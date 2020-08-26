## all my snippets

message(">>> ALL SNIPPETS...")

#### function to enumerate list,  of,  element AND the last one... ----
message("paste_and >>> paste, element and last one")
paste_and <- function(x) {
  return(paste(paste0(x[-length(x)], collapse = ", "), "and", x[length(x)]))
}

message("paste_and_withquote >>> `paste`, `element` and `last one` with quote")
paste_and_withquote <- function(x) {
  return(paste0('`', paste0(x[-length(x)], collapse = "`, `"), "` and `", x[length(x)], "`"))
}

#### head with only nbcol columns show ----
message("h >>> h with nbcol columns")
h <- function(df, nbcol = 5) {
  head(df[, 1:nbcol])
}

#### alias of as.data.frame function----
message("adf >>> alias of as.data.frame")
adf <- as.data.frame

