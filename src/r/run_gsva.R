# Run the pathways (or sets) scores computation
# johaGL 2024
#  -------------
# usage: run_gsva.R DATADIR MATRIX GMT MINSIZE
# example:
# micromamba activate MYMAMBA_OR_CONDA_ENV 
# Rscript zonist/src/r/run_gsva.R $DIR/data/$PRJ/ MX.csv gmt_mz.json 2 

# By defaut this script log transforms the input data: (log(X + 1))

library("GSVA")
library("rjson")


cat("\nThis script runs GSVA (credits: GSVA authors and R package)\n\n")


# functions

run_minimal_test <- function(){
  cat("*-*-*\n DEBUG purpose: running fast example (as no options specified)\n")
  X = matrix(data=rnorm(10*5), nrow=10, ncol=5)
  rownames(X) = as.character(seq(10))
  colnames(X) = as.character(seq(5))
  sets_list = list("set1"=c("5","7","9","6"),
                   "set2"=c("8","9","10"),
                   "set3"=c("1", "10","7"))
  print(sets_list)
  cat("\n printing the fast example (a matrix of scores per set)\n")
  gsvaPar <- gsvaParam(X, sets_list, minSize=2)  
  gsva_scores_matrix <- gsva(gsvaPar, verbose=TRUE)
  print(gsva_scores_matrix)
  cat("\nend of fast example for debug  \n*-*-*\n")
}

# start

args <- commandArgs(trailingOnly = TRUE)     # Note : args[0] # not usable

if (is.na(args[1])){run_minimal_test()}

# variables

data_dir_abspath = args[1]     # DATADIR (absolute/path/to/data)
q_matrix_file = args[2]        # MATRIX  (quantifications)
gmt_file = args[3]             # GMT (sets lists: a.k.a pathways or signatures)
minsize_usr = args[4]          # MINSIZE (minimal nb of features, per pathway)

cat("User options for this analysis:\n")
cat(" data folder (absolute path) : ", data_dir_abspath, "\n")
cat(" quantifications matrix (filename) : ", q_matrix_file, "\n")
cat(" gmt listing the sets (filename) : ", gmt_file, "\n")

json_path <- file.path(data_dir_abspath, gmt_file)
matrix_path <- file.path(data_dir_abspath, q_matrix_file)
out_file = paste0(unlist(strsplit(x=q_matrix_file, split="\\."))[1], 
                  "_gsva.tsv")  

# open files
cat("\nReading files\n")
gmt <- fromJSON(file=json_path)
cat("loaded gmt file with: ", length(gmt), "pathways\n")

cat("loading quantifications...\n")
df <- read.table(matrix_path, row.names=1, header=TRUE, check.names = FALSE,
                sep="\t", comment.char="", 
                stringsAsFactors=FALSE)

X <- as.matrix(df)

colnames(X) = as.character(colnames(X))
cat("loaded quantifications matrix of shape: ", dim(X), "\n")

# compute
cat("Transpose matrix\n")
X <- t(X)
cat("Log transform matrix\n")
X <- log10(X + 1)  # +1 to avoid log zero inf

cat("Computing GSVA ...\n")

gsvaPar <- gsvaParam(X, gmt,
                     minSize=as.integer(minsize_usr))  
rm(X)
gsva_scores_matrix <- gsva(gsvaPar, verbose=FALSE)

#gsva_scores_matrix <- cbind(
#  data.frame("pathway_id"=rownames(gsva_scores_matrix)), 
#  data.frame(gsva_scores_matrix, row.names=NULL))

pathways_ids <- rownames(gsva_scores_matrix) 
print("first pathway:")
print(pathways_ids[1])

exitr <- t(gsva_scores_matrix)

exitr <- cbind(
  data.frame("pathway_id"=rownames(exitr)),
  data.frame(exitr, row.names=NULL)
)
#print(exitr[1:5,1:5])
saveRDS(exitr, 
        file.path(data_dir_abspath, paste0(out_file,".RDS")))

write.table(exitr, 
            file=file.path(data_dir_abspath, out_file),
            sep="\t", row.names=FALSE, col.names=TRUE)

cat("Success ! \n")
cat("Output tab delimited file saved to: ", 
    file.path(data_dir_abspath, out_file) , "\n")
