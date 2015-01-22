# this is the R header to be imported to enable useful packages and self-defined functions.
options(show.error.messages=TRUE)
suppressPackageStartupMessages({  
try(library(MASS))	# standard, no need to install
try(library(cluster))	
try(library(doMC))	# parallele
try(library(stringr)) # string process
try(library(plyr)) # data split and merge
try(library(gdata)) # enforced plot functions
try( library(doSNOW) )	# load doSNOW parallel backend
try( library(lattice) )		# lattice plot package
try( library(geepack) )		# lattice plot package
try( library(data.table) )		# lattice plot package
try( library(xtable) )		# lattice plot package
try( library(getopt) )		# lattice plot package
try( library(dplyr) )		# lattice plot package
try( library(pipeR) )		# lattice plot package
try( library(mvtnorm) )		# lattice plot package
})
options(stringsAsFactors=F)
try(registerDoMC(detectCores()) )
try(options(cores=detectCores()))
# > getOption("cores")
# [1] 8

# > getDoParWorkers()
# [1] 6
#logging_it <- function(...) {
#  catt(...)
#  writeLines(con = )
#}
#---function_that_might_return_null() %||% default value------#
`%||%` <- function(a, b) if (!is.null(a)) a else b

# 0.the two profiling functions. Prefer lineprof.
collect_Rprof <- function(text, output="profile1.out", line.profiling=TRUE){
  Rprof(output, line.profiling = line.profiling)
  eval(parse(text = text, keep.source=TRUE) )	
  Rprof(NULL)
  #---now get summary btw---#
  summaryRprof(output,line="both")
}
collect_lineprof <- function(call){
  l <- lineprof(code = call)
  ll <- align(l)
  arrange(ll,desc(time)) 
}

# 1.
collect_garbage <-function(){while (gc()[2,4] != gc()[2,4] | gc()[1,4] != gc()[1,4]){}}
# 2.
catt <- function(...)
{
 cat(...,'\n' )
}
# 3.
naPlus <- function(x,y) {
    x[is.na(x)] <- 0;
    y[is.na(x)] <- 0;
    x + y 	# this allow pairwise +, not sum to a scalar.
}
# 4.
iblkcol <- function(a, chunks) {
  n <- ncol(a)
  i <- 1

  nextElem <- function() {
    if (chunks <= 0 || n <= 0) stop('StopIteration')
    m <- ceiling(n / chunks)
    r <- seq(i, length=m)
    i <<- i + m
    n <<- n - m
    chunks <<- chunks - 1
    a[,r, drop=FALSE]
  }

  structure(list(nextElem=nextElem), class=c('iblkcol', 'iter'))
}

nextElem.iblkcol <- function(obj) obj$nextElem()
#-------------------#
obj_size <- function(x) print(utils::object.size(x),units='Mb')
save_Rdata <- function(string_obj) save(list=string_obj,file=sprintf('%s.Rdata',string_obj))
read_table <- function(con,...) read.table(con,header=T,check.names=F,sep='\t',fill=T,skip=0,quote='',comment='',strip.white=T,...)
write_matrix <- function(dataframe,filename,...) write.table(dataframe, file = filename,quote=F,row.names=T, col.names=NA,sep='\t',...)
# check history of some pattern for previous cmds.
hishow <- function(pattern) history(pattern=pattern)
# unit is MB here.
purge_bigOBJ <- function(big_cutoff = 100)	# list the objects in memory with corresponding size.
{
    all_obj_strings <- ls(pos = 1)
    # a <- c()
    # for (i in 1:length(all_obj_strings)) {
        # a[[i]] <- cbind(obj = all_obj_strings[i], size = round(object.size(eval(parse(text = all_obj_strings[i])))/1024/1024,
            # 3))
    # }
	all_size = sapply(ls(pos=1),function(x)
                 object.size(get(x)))
                 # format(object.size(get(x)),units="b" ) )
	# format(all_size,units="Mb")
	report = data.frame(names(all_size),all_size,row.names=NULL)
    colnames(report) <- c("obj", "size(mb)")
	report = report[order(as.numeric(report[, 2]),
        decreasing = T), ]
	report[,2] <- round(report[,2]/1024/1024,3)
	return(report)
    # obj_list <- do.call(rbind, a)
    # obj_list_sorted <- obj_list[order(as.numeric(obj_list[, 2]),
        # decreasing = T), ]
    # print(obj_list_sorted, quote = F)
    # invisible(obj_list_sorted)
}

#---------------#
# try to capitalize 1st letter of each word.
simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
      sep="", collapse=" ")
}


nowTimeToString <- function(){
 paste(unlist(str_split(Sys.time(),pattern=' ')),collapse='_')
}

head6 <-function(x,...,ncol=6){
  head(x=x,...)[,1:ncol]
}


head10 <- function(x,...,ncol=10){
    head(x=x,...)[,1:ncol]
  }
 ##================================================================##
 ###  In longer simulations, aka computer experiments,            ###
 ###  you may want to                                             ###
 ###  1) catch all errors and warnings (and continue)             ###
 ###  2) store the error or warning messages                      ###
 ###                                                              ###
 ###  Here's a solution  (see R-help mailing list, Dec 9, 2010):  ###
 ##================================================================##
 tryCatch.W.E <- function(expr)	# for 
{
    W <- NULL
    w.handler <- function(w){ # warning handler
      W <<- w
      invokeRestart("muffleWarning")
    }
    list(value = withCallingHandlers(tryCatch(expr, error = function(e) e),
                                   warning = w.handler),
       warning = W)
}


pchShow <-	# demonstrate all the symbols for plot convenience.
       function(extras = c("*",".", "o","O","0","+","-","|","%","#"),
                cex = 3, ## good for both .Device=="postscript" and "x11"
                col = "red3", bg = "gold", coltext = "brown", cextext = 1.2,
                main = paste("plot symbols :  points (...  pch = *, cex =",
                             cex,")"))
       {
         nex <- length(extras)
         np  <- 26 + nex
         ipch <- 0:(np-1)
         k <- floor(sqrt(np))
         dd <- c(-1,1)/2
         rx <- dd + range(ix <- ipch %/% k)
         ry <- dd + range(iy <- 3 + (k-1)- ipch %% k)
         pch <- as.list(ipch) # list with integers & strings
         if(nex > 0) pch[26+ 1:nex] <- as.list(extras)
         plot(rx, ry, type="n", axes = FALSE, xlab = "", ylab = "",
              main = main)
         abline(v = ix, h = iy, col = "lightgray", lty = "dotted")
         for(i in 1:np) {
           pc <- pch[[i]]
           ## 'col' symbols with a 'bg'-colored interior (where available) :
           points(ix[i], iy[i], pch = pc, col = col, bg = bg, cex = cex)
           if(cextext > 0)
               text(ix[i] - 0.3, iy[i], pc, col = coltext, cex = cextext)
         }
       }

# > dim(datExp_list_numeric_only8cancer_transposedGeneAsColumns$Breast_datExp)
# [1]   529 17813
# try(allowWGCNAThreads(8))
# Allowing multi-threading with up to 8 threads.

# for loop
#===========Logical predicates and boolean algebra on function level=================#
and <- function(f1, f2) {
force(f1); force(f2)
function(...) {
f1(...) && f2(...)
}
}
or <- function(f1, f2) {
force(f1); force(f2)
function(...) {
f1(...) || f2(...)
}
}
not <- function(f1) {
force(f1); force(f2)
function(...) {
!f1(...)
}
}
# example:
# Filter(or(is.character, is.factor), iris)
options(show.error.messages=TRUE)
