###optmatch:::fullmatchNumControlsHandling
try(optmatch:::fullmatchNumControlsHandling(letters, letters, whicharg="min.controls"))
optmatch:::fullmatchNumControlsHandling(1:2, character(0), whicharg="min.controls")
optmatch:::fullmatchNumControlsHandling(1:2, "m", whicharg="min.controls")
try(optmatch:::fullmatchNumControlsHandling(1:2, letters, whicharg="min.controls"))
try(optmatch:::fullmatchNumControlsHandling(c(A=1, B=2), letters, whicharg="min.controls"))
optmatch:::fullmatchNumControlsHandling(2, letters[1:3], whicharg="min.controls")
optmatch:::fullmatchNumControlsHandling(c(a=1,b=2,c=3), letters[1:3], whicharg="min.controls")
optmatch:::fullmatchNumControlsHandling(with(expand.grid(A=1:2, B=letters[1:3]), tapply(A, B, mean)), letters[1:3], whicharg="min.controls")



