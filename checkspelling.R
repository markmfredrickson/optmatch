check_spelling <- function(path) {
  # this is usually loaded anyway, but just in case...
  library("utils")

  words <- aspell_package_Rd_files(path, program = "aspell", control = list("-p ./lexicon.txt"))

  if (dim(words)[1] > 0) {
    printAspellByFile(words)
    stop("Correct spelling mistakes before proceeding.")
  }

  message("No spelling errors detected in documentation.")
}

printAspellByFile <- function (x, sort = TRUE, verbose = FALSE, indent = 2L, ...) 
{
    if (!(nr <- nrow(x))) 
        return(invisible(x))
    if (sort) 
        x <- x[order(x$File, x$Original, x$Line, x$Column), ]
    if (verbose) 
        out <- sprintf("%sWord: %s (%s:%d:%d)\n%s", c("", rep.int("\n", 
            nr - 1L)), x$Original, x$File, x$Line, x$Column, 
            formatDL(rep.int("Suggestions", nr), sapply(x$Suggestions, 
                paste, collapse = " "), style = "list"))
    else {
        s <- split(sprintf("%s:%d:%d", x$Original, x$Line, x$Column), 
            x$File)
        sep <- sprintf("\n%s", paste(rep.int(" ", indent), collapse = ""))
        out <- paste(names(s), sapply(s, paste, collapse = sep), 
            sep = sep, collapse = "\n\n")
    }
    writeLines(out)
    invisible(x)
}

make_dictionary <- function(path) {
  words <- aspell_package_Rd_files(path, program = "aspell")
  aspell_write_personal_dictionary_file(words, out = "lexicon.txt", program = "aspell")  
}
