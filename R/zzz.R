.onAttach <- function(lib, pkg) {   
packageStartupMessage(paste(
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
"fullmatch() function. For more information, enter\n",
"relaxinfo() at the command line.\n"
))
} 
