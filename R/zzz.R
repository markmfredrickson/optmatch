
.onLoad <- function(lib, pkg) {   
if (!interactive())
  {
cat(
"You're loading optmatch, by Ben Hansen, a package for flexible\n",
"and optimal matching.  Important license information:\n",
"The optmatch package makes essential use of D. P. Bertsekas\n",
"and P. Tseng\'s RELAX-IV algorithm and code, as well as\n" ,
"Bertsekas\' AUCTION algorithm and code.\n",
"Bertsekas and Tseng freely permit their software to be used for\n",
"research purposes, but non-research uses, including the use of it\n",
"to \'satisfy in any part commercial delivery requirements to\n",
"government or industry,\' require a special agreement with them.\n",
"By extension, this requirement applies to any use of the \n",
"fullmatch() function. (If you are using another package that has\n",
"loaded optmatch, then you will probably be using fullmatch indirectly.)\n",
"For more information, enter relaxinfo() at the command line\n"
)
} else
{
outFile <- tempfile()
outConn <- file(outFile, open="w")
writeLines(paste(
"You're loading optmatch, a package for flexible, optimal matching.\n\n ",
"Important license information: This package makes use of Bertsekas\nand ",
"Tseng\'s RELAX-IV algorithm and code, as well as Bertsekas\'\nAUCTION " ,
"algorithm and code.  Bertsekas and Tseng freely permit\ntheir ",
"software to be used for research purposes, but non-research\nuses, ",
"including uses to \'satisfy in any part commercial delivery\nrequirements ",
"to government or industry,\' require a special agreement.\nBy ",
"extension, this requirement applies to most any use of the \nfullmatch() ",
"function. (If you are using another package that has\nloaded optmatch, ",
"then you may be using fullmatch indirectly.)\n\nFor ",
"more information, enter relaxinfo() at the command line.\n To ",
"dismiss this message, press \'q\'. -B. Hansen",
sep=""), 
outConn)
close(outConn)
file.show(outFile, 
title="Important license information regarding optmatch package",
delete.file=TRUE) 
}
}
