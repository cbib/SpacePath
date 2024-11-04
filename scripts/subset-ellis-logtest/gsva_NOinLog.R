# Run the pathways (or sets) scores computation
# johaGL 2024
#  -------------
# Assuming you have already created micromamba 'r_gsva' environment:
# micromamba activate r_gsva

# ** ACCEPTS TRANSPOSED (tab delim) : m_pixels (rows) x n_features (columns) **

#  (not using logarithm for the input)

library("GSVA")
library("rjson")


cat("\nThis script runs GSVA (credits: GSVA authors and R package)\n\n")

# functions

run_minimal_test <- function(){
  cat("*-*-*\n DEBUG purpose: running fast example\n")
  X = matrix(data=rnorm(10*5), nrow=10, ncol=5)
  rownames(X) = as.character(seq(10))
  colnames(X) = as.character(seq(5))
  sets_list = list("set1"=c("5","7","9","6"),
                   "set2"=c("8","9","10"),
                   "set3"=c("1", "10","7"))
  print(sets_list)
  cat("\n printing the fast example (a matrix of scores per set)\n")
  gsvaPar <- gsvaParam(X, sets_list)  # defaults
  gsva_scores_matrix <- gsva(gsvaPar, verbose=FALSE)
  print(gsva_scores_matrix)
  cat("\nend of fast example for debug  \n*-*-*\n")
}

# start
args <- commandArgs(trailingOnly = TRUE)     # Note : args[0] # not usable
# run_minimal_test()   # un-comment only if debug required

# variables

data_dir_abspath = args[1]     # DATADIR (absolute/path/to/data)
q_matrix_file = args[2]        # MATRIX  (quantifications)
gmt_file = args[3]             # GMT (sets lists: a.k.a pathways or signatures)

cat("User options for this analysis:\n")
cat(" data folder (absolute path) : ", data_dir_abspath, "\n")
cat(" quantifications matrix (filename) : ", q_matrix_file, "\n")
cat(" gmt listing the sets (filename) : ", gmt_file, "\n")

json_path <- file.path(data_dir_abspath, gmt_file)
matrix_path <- file.path(data_dir_abspath, q_matrix_file)
out_file = paste0(unlist(strsplit(x=q_matrix_file, split="\\."))[1], 
                  "_gsva-NOinLOG.tsv")   #  output file name is here automatic
# TODO: change the output name

# open files
cat("\nReading files\n")
gmt <- fromJSON(file=json_path)
cat("loaded gmt file with: ", length(gmt), "pathways\n")

cat("loading quantifications...\n")

df <- read.table(matrix_path, row.names=1, header=TRUE, 
                 check.names = FALSE,
                 sep="\t", comment.char="", 
                 stringsAsFactors=FALSE)

X <- as.matrix(df)
rm(df)
colnames(X) = as.character(colnames(X))
cat("loaded quantifications matrix of shape: ", dim(X), "\n")


# compute
cat("Transpose matrix\n")
X <- t(X)

cat("Computing GSVA ...\n")

gsvaPar <- gsvaParam(X, gmt)  # run with defaults
rm(X)
gsva_scores_matrix <- gsva(gsvaPar, verbose=FALSE)

pathways_ids <- rownames(gsva_scores_matrix) 
print("first pathway:")
print(pathways_ids[1])

exitr <- t(gsva_scores_matrix)

exitr <- cbind(
  data.frame("pathway_id"=rownames(exitr)),
  data.frame(exitr, row.names=NULL)
)
#print(exitr[1:5,1:5])

write.table(exitr, 
            file=file.path(data_dir_abspath, out_file),
            sep="\t", row.names=FALSE, col.names=TRUE)

cat("Success ! \n")
cat("Output tab delimited file saved to: ", 
    file.path(data_dir_abspath, out_file) , "\n")
