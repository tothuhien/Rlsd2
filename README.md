# Rlsd2
This package is developped as an R-wrapper for lsd2 (https://github.com/tothuhien/lsd2).
For now, the package is only tested under Unix system.

All options and file format used here are the same as the original lsd2, but this R package also accepts some R-object inputs. Users could refer to that the original lsd2 for better documentation.

To install it (some packages are required: devtools, ape, treeio): 

    devtools::install_github("tothuhien/Rlsd2")

Example of using it:

    library(Rlsd2)
    result <- lsd2(inputTree="data/D750_11_10_rooted.tree", 
                   inputDate="data/D750_11_10.date", 
                   outFile = "data/test_lsd2", seqLen=1000)
                 
Another example with R-object input

    tree <- read.tree("data/D750_11_10_rooted.tree")
    dateTbl <- read.table("data/D750_11_10.date",skip=1,colClasses = "character")
    date <- dateTbl[,2]
    names(date) <- dateTbl[,1]
    result <- lsd2(inputTree=tree, inputDate=date, 
                   outFile = "data/test_lsd2", seqLen=1000)

`?lsd2` for help page

