##' Function to display license information regarding code embedded in
##' \code{optmatch}.
##'
##' @title Display license information about embedded code
##' @return Some information about licenses of code and algorithms on
##'   which \code{fullmatch} depends.
##' @author Ben B. Hansen
##' @export
relaxinfo <- function() {
  if (!interactive()) {
    cat(
      "The optmatch package makes essential use of D. P. Bertsekas\n",
      "and P. Tseng\'s RELAX-IV algorithm and code, as well as\n" ,
      "Bertsekas\' AUCTION algorithm and code:\n\n",
      "Bersekas, D. P. and  Tseng, P., \"Relaxation Methods for\n",
      "Minimum Cost ...\" Operations Research, vol. 26, 1988, 93-114;\n",
      "Bertsekas, D. P., \"An Auction/Sequential Shortest Path\n" ,
      "Algorithm for the Minimum Cost Flow Problem\", LIDS Report\n" ,
      "P-2146, MIT, Nov. 1992;\n",
      "<http://web.mit.edu/dimitrib/www/noc.htm>.\n\n",
      "Bertsekas and Tseng freely permit their software to be used for\n",
      "research purposes, but non-research uses, including the use of it\n",
      "to \'satisfy in any part commercial delivery requirements to\n",
      "government or industry,\' require a special agreement with them.\n",
      "By extension, this requirement applies to most any use of R\n",
      "functions in the optmatch package.\n\n",
      "To request permissions not here relayed, contact Professor Bertsekas at\n",
      "Laboratory for Information and Decision Systems\n",
      "Massachusetts Institute of Technology\n",
      "Cambridge, MA 02139\n",
      "(617) 253-7267 <dimitrib@mit.edu>\n"
    )
  } else {
    outFile <- tempfile()
    outConn <- file(outFile, open="w")
    writeLines(paste(
      "The optmatch package makes essential use of D. P. Bertsekas\nand ",
      "P. Tseng\'s RELAX-IV algorithm and code, as well as\nBertsekas\' " ,
      "AUCTION algorithm and code:\n\nBersekas, ",
      "D. P. and  Tseng, P., \"Relaxation Methods for\nMinimum ",
      "Cost ...\" Operations Research, vol. 26, 1988, 93-114;\nBertsekas, ",
      "D. P., \"An Auction/Sequential Shortest Path\nAlgorithm " ,
      "for the Minimum Cost Flow Problem\", LIDS Report\nP-2146, " ,
      "MIT, Nov. 1992;\n<http:",
      "//web.mit.edu/dimitrib/www/noc.htm>.\n\nBertsekas ",
      "and Tseng freely permit their software to be used for\nresearch ",
      "purposes, but non-research uses, including the use of it\nto ",
      "\'satisfy in any part commercial delivery requirements to\ngovernment ",
      "or industry,\' require a special agreement with them.\nBy ",
      "extension, this requirement applies to most any use of R\nfunctions ",
      "in the optmatch package.\n\nTo ",
      "request permissions not here relayed, contact Professor ",
      "Bertsekas at\nLaboratory ",
      "for Information and Decision Systems\nMassachusetts ",
      "Institute of Technology\nCambridge, ",
      "MA 02139\n(617) 253-7267 <dimitrib@mit.edu>\n\n(To ",
      "dismiss this message, press \'q\'.)",
      sep=""),
      outConn)
    close(outConn)
    file.show(outFile,
              title="Important license information regarding optmatch package",
              delete.file=TRUE)
  }
}
