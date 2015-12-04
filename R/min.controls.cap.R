#' @export
minControlsCap <- function(distance, max.controls=NULL)
{
  distance <- as.matrix(distance) # cast ISM to matrix, temporary
  if (!is.list(distance) & !is.matrix(distance))
    stop("argument \'distance\' must be a matrix or list")

  if (is.matrix(distance))
  {
    tmp <- maxControlsCap(t(distance),  min.controls=
                                          switch(1+is.null(max.controls),
                                                 ifelse(max.controls>=1, 1/ceiling(max.controls),
                                                        floor(1/max.controls) ), NULL))
  } else   {
    tmp <- maxControlsCap(lapply(distance, t), min.controls=
                                                 switch(1+is.null(max.controls),
               ifelse(max.controls>=1, 1/ceiling(max.controls),
                      floor(1/max.controls) ), NULL))
  }

  list(strictest.feasible.min.controls=
         (1/tmp$strictest.feasible.max.controls),
       given.max.controls=(1/tmp$given.min.controls))
}
