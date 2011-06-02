fullmatchNumControlsHandling <- function(numcontrols,idclevels,whicharg)
  {

    stopifnot(whicharg %in% c('min.controls', 'max.controls'))
    
    if (is.null(numcontrols))
      numcontrols <- c(0,Inf)[match(whicharg, c('min.controls', 'max.controls'))]

lennumcontrols <- length(numcontrols)

if (!is.numeric(numcontrols))
  stop(paste("is.numeric(",whicharg ,")", " is not TRUE"))
if (any(is.na(numcontrols)))
  stop(paste("NAs or NaNs in argument",whicharg))

if (lennumcontrols>1 && length(idclevels)<=1)
  {
    numcontrols <- numcontrols[1]
    lennumcontrols <- 1
    warning(paste("Ignoring elements of", whicharg, "other than first, b/c I see no subclasses."))
  }
if (lennumcontrols>1 && is.null(names(numcontrols)))
  stop(paste("Argument", whicharg, "needs names to align with subclasses."))

if (lennumcontrols>1 && !all(idclevels %in% names(numcontrols)))
  stop(paste("'",whicharg ,"'"," not specified for all subclasses.", sep=""))

if (lennumcontrols==1)
  {
   ans <- rep(numcontrols, max(1,length(idclevels)))
   if (length(idclevels)) names(ans) <- idclevels
  } else ans <- numcontrols[idclevels]
ans
  }
