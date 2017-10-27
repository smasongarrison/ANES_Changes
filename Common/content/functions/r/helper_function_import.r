library(RCurl)
setwd(tempdir())

destfile = "MasonHelperFunctions.R"
x = getBinaryURL("https://dl.dropbox.com/s/5tsqgothgdfbage/MasonHelperFunctions.R?dl=0", followlocation = TRUE, ssl.verifypeer = FALSE)
writeBin(x, destfile, useBytes = TRUE)
source(paste(tempdir(), "/MasonHelperFunctions.R", sep = ""))

# remove files from tempdir:
unlink(dir())