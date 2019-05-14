setClass("MatchablesInfo", contains="data.frame")
setValidity("MatchablesInfo", function(object){
    errors <- character(0)
    for (nm in c('name', 'subproblem', 'kind'))
    {
        if ( !(nm %in% colnames(object)) )
            errors  <- c(errors,
                         paste0("Need a '", nm, "' column.")
                         )
        if ( !(is.character(object[[nm]])) )
            errors  <- c(errors,
                         paste0("Column '", nm, "' should have type character.")
                         )
    }
    if ( !all(object[['kind']] %in% c("treatment", "control")) )
        errors  <- c(errors,
                     "'kind' values other than 'treatment' or 'control'.")
    if (length(errors)==0) TRUE else errors  
})
