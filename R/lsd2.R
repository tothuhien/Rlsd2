#' @title A phylogeny dating method using least-squares algorithms and criteria
#' @description This an R-wrapper of the tool lsd2 (https://github.com/tothuhien/lsd2). Users could refer to that tool for better documentation.
#' All options and file format used here are the same as the original lsd2, but this R package also accepts some R-object inputs.
#' @param inputTree either the input file containing tree(s) in newick format, or a phylo or multiphylo tree object.
#' @param inputDate either the input date file, or a vector of date whose names correpond to the names of the dated node.
#' @param seqLen the length of sequences that were used to build the tree(s)
#' @param partitionFile the input partition file if there's any
#' @param outFile the basename of the output files, if NA then lsd2 will generate temporary files within working folder.
#' @param outGroup either the file contains the outgroups, or a vector of outgroup names. The tree will be rooted based on the specified outgroups. Note that the outgroups must form a monophyletic in the input tree.
#' @param givenRate either the file contains the given rates, or a vector of given rates, each rate corresponds to a tree in the input tree(s)
#' @param constraint if TRUE then temporal constraints is imposed
#' @param variance apply variances for branche lengths. Either 0, or 1, or 2 which respectively means without variance, run lsd2 once with variances, and run lsd2 twice with variances where the second time variances are based on the estimated branch length of the first run.
#' @param confidenceInterval if specified, then will compute the confidence intervals. 
#' This can be either the number of bootstraps to generate simulated trees for calculing confidence intervals or a file contains the bootstrap trees or a list of bootstrap trees. 
#' Note that those bootstrap trees have to have exactly the same topology as the input tree(s).
#' The confidence interval date of each node is reported in the nexus outfile.
#' @param splitInternalExternalBranches whether to estimate the rate of internal and external branches separately. False by default.
#' @param minRate the minimum threshold for rate
#' @param estimateRoot either "l","a" or "as", to esimate or re-estimate the roots. Use "l" (only for rooted tree) to re-estimate the root around the given root, "a" to (re)-estimate the root over all branches using fast methode, "as" to (re)-esimate to root over all branches with slow methode (apply temporal constraints through-out the whole search process)
#' @param b the parameter to adjust variance, estimate by default
#' @param q the standard deviation of lognormal distribution to apply on branch lengths when calculating confidence intervals. 0.2 by default.
#' @param rootDate the date of the root. 
#' @param tipsDate the date of all leaves
#' @param verbose verbose mode for output messages.
#' @param keepOutgroup TRUE to keep the outgroups specifed in `outGroup`, FALSE to remove them.
#' @param nullblen every branch length <= `nullblen` will be collapsed before dating. By default is 0.5/`seqLen`
#' @param support if the tree containts support value and you want to collapse all branches having small support value, then you can specify the thresold of support here.
#' @param minblen the minimum length for internal branches of the output time-scaled tree. 
#' Can be either a positive real, or the letter "e" for letting the program to estimate this value.
#' minblen is rounded to a meaningful unit (for example date, week, year ...) using `roundTime` parameter.
#' 0 by default.
#' @param minblenL the minimum length for external branches of the output time-scaled tree. Similar to minblen.
#' @param roundTime the factor to round the time of minblen and minblenL. By default, factor 365 rounds the minblen to number of days.
#' @param outDateFormat the format of the output date, 1 for real, and 2 for year-month-day. By default, the output date format is based on the input date format.
#' @param m the number of sampling dated nodes to calculate median rates, used in estimating outliers nodes and minblen. By default is 10.
#' @param ZscoreOutlier if specify then lsd2 will esimate and remove outliers nodes before dating. A normal value of ZscoreOutlier could be 3, but you can adjust it bigger/smaller depending if you want to have less/more outliers. Note that for now, some functionalities could not be combined with outliers estimation, for example estimating multiple rates, imprecise date constraints.
#' @param nbData the number of tree(s) in the `inputTree` that are going to process
#' @examples
#' result <- lsd2(inputTree="data/D750_11_10_rooted.tree", inputDate="data/D750_11_10.date", outFile = "data/test_lsd2", seqLen=1000)
#' ## or
#' tree <- read.tree("data/D750_11_10_rooted.tree")
#' dateTbl <- read.table("data/D750_11_10.date",skip=1,colClasses = "character")
#' date <- dateTbl[,2]
#' names(date) <- dateTbl[,1]
#' result <- lsd2(inputTree=tree, inputDate=date, outFile = "data/test_lsd2", seqLen=1000)
#' @importFrom ape read.tree write.tree
#' @importFrom treeio read.beast
#' @return a list that includes: estimated rate(s), root date(s), estimated tree(s), and the output file names.
#' @export
#'
#'
lsd2 <- function(inputTree, inputDate = NA, seqLen, partitionFile = NA, outFile = NA, outGroup = NA, givenRate = NA,
                 constraint = TRUE, variance = 1, confidenceInterval = NA, splitInternalExternalBranches = FALSE,
                 minRate = 1e-10, estimateRoot = NA, b = NA, q = 0.2,
                 rootDate=NA, tipsDate =NA, verbose = FALSE, keepOutgroup = TRUE,
                 nullblen = NA, support =  NA,  minblen = NA, minblenL = NA,
                 roundTime = NA, outDateFormat = NA,  m = 10, ZscoreOutlier = NA, nbData = 1){
  if ((typeof(inputTree) == "character") && file.exists(inputTree)){
    inputTree = normalizePath(inputTree)
    if (is.na(outFile)) outFile = paste0(inputTree,".result")
  } else{
    f = tempfile()
    if (class(inputTree)=="phylo" || class(inputTree)=="multiPhylo"){
      write.tree(inputTree,f)
      inputTree=f
      if (is.na(outFile)) outFile = tempfile(pattern = "lsd2_",tmpdir = getwd());
    } else {
      cat("inputTree is not recognized as an exist file or phylo/multiPhylo object\n")
      return(NULL)
    }
  }
  if (!is.na(inputDate)){
    if ((typeof(inputDate) == "character") && file.exists(inputDate)){
      inputDate = normalizePath(inputDate)
    } else {
      d = tempfile()
      cat(length(inputDate),"\n",file = d)
      for (i in 1:length(inputDate)){
        cat(names(inputDate)[i],inputDate[i],"\n",append = T,file = d)
      }
      inputDate = d
    }
  }
  if (!is.na(outGroup)){
    if ((typeof(outGroup) == "character") && file.exists(outGroup)){
      outGroup = normalizePath(outGroup)
    } else {
      outgroup = tempfile()
      if (!is.na(outGroup)){
        cat(length(outGroup),"\n",file = outGroup)
        for (i in 1:length(outGroup)){
          cat(outGroup[i],"\n",append = T,file = outgroup)
        }
        outGroup = outgroup
      }
    }
  }
  if (!is.na(givenRate)){
      if ((typeof(givenRate) == "character") && file.exists(givenRate)){
        givenRate = normalizePath(givenRate)
      } else {
        rateFile = tempfile()
        if (!is.na(givenRate)  && typeof(givenRate)=="double"){
          for (i in 1:length(givenRate)){
            cat(givenRate[i],"\n",append = T,file = rateFile)
          }
          givenRate = rateFile
        }
      }
  }
  if (!is.na(partitionFile)){
    if ((typeof(partitionFile) == "character") && file.exists(partitionFile)){
      partitionFile = normalizePath(partitionFile)
    }
  }
  if (!is.na(confidenceInterval)){ 
    if (class(confidenceInterval)=="multiPhylo"){
      f = tempfile()
      write.trees(confidenceInterval,f)
      confidenceInterval=f
    } else if (class(confidenceInterval) != "character"){
      confidenceInterval=as.integer(confidenceInterval)  
      if (is.na(confidenceInterval)){
        cat("confidenceInterval must be either a positive integer, or a file containing bootstrap trees, or a multiPhylo object\n")
        return(NULL)
      }
    } else{
      confidenceInterval = normalizePath(confidenceInterval)
      if (!file.exists(confidenceInterval)){
        cat("The bootstrap trees file does not exist.\n")
        return(NULL)
      }
    }
  }
  if (!is.logical(splitInternalExternalBranches)){
    cat("splitInternalExternalBranches must be either TRUE or FALSE\n")
    return(NULL)
  }
  if (!is.na(outDateFormat)) outDateFormat=as.integer(outDateFormat)
  if (!is.numeric(variance) || (variance!=0 && variance!= 1 && variance != 2)) {
    cat("Variance must be either 0, 1, or 2\n")
    return(NULL)
  }
  if (!is.na(m)) m=as.integer(m)
  if (!is.na(nbData)) nbData=as.integer(nbData)

  res = .Call("Rlsd2",inputTree,  inputDate,  partitionFile,  outFile,  outGroup,  givenRate,
        seqLen,  constraint ,  variance ,  confidenceInterval , splitInternalExternalBranches, minRate ,
        estimateRoot ,  b ,  q ,  rootDate ,  tipsDate ,  verbose ,  keepOutgroup ,
        nullblen ,  support ,   minblen ,  minblenL ,
        roundTime,  outDateFormat , m , ZscoreOutlier, nbData)


  if (!is.null(res)){
    names(res) <- c("rate","tMRCA")
    cat("Reading output trees ...\n")
    res[["newickTree"]] = read.tree(paste0(outFile,".nwk"))
    res[["nexusTreeFile"]] = read.beast(paste0(outFile,".nexus"))
    res[["dateNexusTreeFile"]] = read.beast(paste0(outFile,".date.nexus"))
    res[["outResultFiles"]] = c(outFile,paste0(outFile,".nwk"),paste0(outFile,".nexus"),paste0(outFile,".date.nexus"))
    cat("Done.\n")
  }
  return (res)
}
