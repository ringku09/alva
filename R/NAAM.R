library(stringr)

# prepare corpus
corpus <- matrix(c(1,2,0,0,0,0,
                   3,1,0,0,0,0,
                   2,0,0,0,0,0,
                   3,3,2,3,2,4,
                   0,0,3,2,0,0,
                   0,0,4,1,0,0,
                   0,0,0,0,4,3,
                   0,0,0,0,2,1,
                   0,0,0,0,3,2,
                   0,0,1,0,2,3), ncol=6, byrow=T)
vv <- naov_data2$compound[gene2, ]
rownames(vv) <- get_genename(rownames(vv), organism = "rat")
vv_mat <-exp(vv)
#GCmat <-round(1/(1+exp(-vv))*100) #vv_mat/sum(vv_mat)
#pp <- PHVM(GCmat,k = 3)
corpus <- vv_mat/sum(vv_mat)  #round(1/(1+exp(-vv))*100)
# corpus <- apply(vv_mat, 2, function(col) col / sum(col))
sum(corpus)
# initialize parameters
ntopic <- 3
#docnames <- c('Doc1','Doc2','Doc3','Doc4','Doc5','Doc6')
#termnames <- c('Baseball','Basketball','Boxing','Money','Interest','Rate','Democrat','Republican','Cocus','President')
docnames <- colnames(corpus)
termnames <- rownames(corpus)
#colnames(corpus) <- docnames
ndocs <- length(docnames)
#rownames(corpus) <- termnames
nterms <- length(termnames)
topicnames <- paste0('topic',1:ntopic)
dtnames <- c()
for (i in 1:dim(corpus)[2]) {
  dtnames <- c(dtnames,paste0(docnames[i],' ',termnames))
}
posterior.init <- matrix(runif(dim(corpus)[1]*dim(corpus)[2]*ntopic,min=0,max=1),ncol=ntopic)
colnames(posterior.init) <- topicnames
rownames(posterior.init) <- dtnames
pz.init <- matrix(runif(ntopic,min=0,max=1),ncol=ntopic)
colnames(pz.init) <- topicnames
pdz.init <- matrix(runif(dim(corpus)[2]*ntopic,min=0,max=1),ncol=ntopic)
colnames(pdz.init) <- topicnames
rownames(pdz.init) <- docnames
pwz.init <- matrix(runif(dim(corpus)[1]*ntopic,min=0,max=1),ncol=ntopic)
colnames(pwz.init) <- topicnames
rownames(pwz.init) <- termnames
parameter.init <- list(pwz.init,pdz.init,pz.init)

# Expectation Step
estep <- function(parameter,posterior) {
  pwz <- parameter[[1]]
  pdz <- parameter[[2]]
  pz <- parameter[[3]]
  for (i in 1:(dim(corpus)[1]*dim(corpus)[2])) {
    doc <- unlist(strsplit(dtnames[i],' '))[1]
    term <- unlist(strsplit(dtnames[i],' '))[2]
    denominator <- sum(pz * pwz[which(rownames(pwz)==term),] * pdz[which(rownames(pdz)==doc),])
    for (j in 1:ntopic) {
      numerator <- pz[1,j] * pdz[which(rownames(pdz)==doc),j] * pwz[which(rownames(pwz)==term),j]
      posterior[i,j] <- numerator/denominator
    }
  }
  return(posterior)
}

# Maximization Step
mstep <- function(posterior, parameter) {
  pwz <- parameter[[1]]
  pdz <- parameter[[2]]
  pz <- parameter[[3]]
  for (i in 1:dim(pwz)[1]) {
    for (j in 1:dim(pwz)[2]) {
      pwznumerator <- sum(corpus[i,] * posterior[which(str_detect(rownames(posterior), termnames[i])),j])
      pwzdenominator <- sum(corpus * posterior[,j])
      pwz[i,j] <- pwznumerator/pwzdenominator
    }
  }

  for (i in 1:dim(pdz)[1]) {
    for (j in 1:dim(pdz)[2]) {
      pdznumerator <- sum(corpus[,i] * posterior[which(str_detect(rownames(posterior), docnames[i])),j])
      pdzdenominator <- sum(corpus * posterior[,j])
      pdz[i,j] <- pdznumerator/pdzdenominator
    }
  }

  for (i in 1:dim(pz)[2]) {
    pznumerator <- sum(posterior[,i] * corpus)
    pzdenominator <- sum(corpus)
    pz[1,i] <- pznumerator/pzdenominator
  }

  return(list(pwz,pdz,pz))
}

# calculate probs
posterior.iter <- estep(parameter.init, posterior.init)
parameter.iter <- mstep(posterior.init, parameter.init)

while(i<100) {
  posterior.iter <- estep(parameter.iter, posterior.iter)
  parameter.iter <- mstep(posterior.iter, parameter.iter)
  i <- i + 1
}

dim(vv)

image(posterior.iter)

mat <- round(1/(1+exp(-abs(vv)))*100)

svs_pls <- fast_psa(corpus, k = 3)

pheatmap(jp, scale = "column", color = hcl.colors(1000, "PiYG"))
image(jp)

jp <- svs_pls$prob1 %*% diag(svs_pls$prob0) %*% t(svs_pls$prob2)

colnames(posterior.iter)[apply(posterior.iter, 1, which.max)]

Y <- jp
for (i in 1:nrow(jp)) {
  sorted_indices <- order(jp[i, ])
  Y[i, ] <- jp[i, sorted_indices]
}

jp <- parameter.iter[[1]] %*% diag(parameter.iter[[3]][1,]) %*% t(parameter.iter[[2]])

colSums(parameter.iter[[2]])
image(parameter.iter[[2]])
pheatmap(posterior.iter, scale = "row", color = hcl.colors(1000, "PiYG"))




#########################     vectorized code  ########################
# Modifies the code in
#    github thread: https://gist.github.com/ratsgo/c25deb6d79f3ab8b0b050af751fbbdb8#file-plsa-r-L78
# to speed it up by introducing vectorisation alternatives to nested for loops
# and by avoiding cumbersome string lookups

# assuming this is based on Hofmann (2001) https://doi.org/10.1023/A:1007617005950, p. 182 eqns 6, 11 and 12 but re-casted for parametrization as in eq. 5


# prepare corpus
corpus <- matrix(c(1,2,0,0,0,0,
                   3,1,0,0,0,0,
                   2,0,0,0,0,0,
                   3,3,2,3,2,4,
                   0,0,3,2,0,0,
                   0,0,4,1,0,0,
                   0,0,0,0,4,3,
                   0,0,0,0,2,1,
                   0,0,0,0,3,2,
                   0,0,1,0,2,3), ncol=6, byrow=T)

# set labels
docnames <- c('Doc1','Doc2','Doc3','Doc4','Doc5','Doc6')
termnames <- c('Baseball','Basketball','Boxing','Money','Interest','Rate','Democrat','Republican','Cocus','President')
colnames(corpus) <- docnames
rownames(corpus) <- termnames

# set counters
nterms <- length(termnames)
ndocs <- length(docnames)
n_docs <- ndocs

# initialize parameters
ntopic <- 3
topicnames <- paste0('topic',1:ntopic)


#### vectorisation setup ####
sum_corpus <- sum(corpus)                               # a constant for later use

# vectorize row and col indexes
a_plsa <- 1:ncol(corpus)
b_plsa <- 1:nrow(corpus)
for_idx_plsa <- expand.grid(a_plsa,b_plsa)              # expand grid pre-generates combinations of indexes to lookup.  in the original thread this step is a nested for loop coupled with cumbersome unlisting
colnames(for_idx_plsa) <- c("doc","semant_unit")

# vectorize again: fix doc, change semantic unit
for_idx_plsa_swapped <- expand.grid(b_plsa,a_plsa)
colnames(for_idx_plsa_swapped) <- c("semant_unit","doc")

# vectorize corpus PLSA to make  it commensurate to the expand.grid
vect_corpus_PLSA <- c(t(corpus))                       # to compute pwz (docs * terms)
vect_corpus_PLSA_2 <- c(corpus)                        # to compute pdz (terms * docs)


#### random initialisation ####
set.seed(123)

pz.init <- matrix(runif(ntopic,min=0,max=1),ncol=ntopic)
colnames(pz.init) <- topicnames

pdz.init <- matrix(runif(ndocs*ntopic,min=0,max=1),ncol=ntopic)
colnames(pdz.init) <- topicnames
rownames(pdz.init) <- docnames

pwz.init <- matrix(runif(nterms*ntopic,min=0,max=1),ncol=ntopic)
colnames(pwz.init) <- topicnames
rownames(pwz.init) <- termnames

parameter.init <- list(pwz.init,pdz.init,pz.init)


#### Expectation Step - without for loops ####
estep2 <-  function(x, y, parameter){

  pwz <- parameter[[1]]
  pdz <- parameter[[2]]
  pz <- parameter[[3]]

  a <- pwz[y,]
  b <- pdz[x,]

  # vectorize pz
  c <- matrix(pz, nrow=nrow(y), ncol=length(pz), byrow=TRUE)         # thread: https://stackoverflow.com/questions/14927507/duplicate-vector-into-matrix-r

  numerator <- c * b * a
  denominator <- apply(numerator,1,sum)

  numerator/denominator
}


#### Maximization Step - without for loops ####
mstep_pwz <-  function(x, posterior){
  temp_result_1 <- apply(posterior,2,function(w){       # posterior also has (docs*terms) rows but the n. of columns equals the n of topics. Each col has to be multiplied
    w*vect_corpus_PLSA                                  # element-wise multiplication. To be summed
  })
  temp_result_2 <- cbind(x, temp_result_1)
  colnames(temp_result_2) <- c("semant_unit", colnames(temp_result_1))
  temp_result_3 <- aggregate(. ~semant_unit, data=as.data.frame(temp_result_2), sum)               # kind of sumif thread: https://stackoverflow.com/questions/21607464/what-is-the-equivalent-of-the-sumif-function-in-r
  pwznumerator <- temp_result_3[,-1]                                                               # remove first column
  pwzdenominator <- apply(pwznumerator, 2, sum)
  apply(pwznumerator,1, function(x){
    x * 1/pwzdenominator
  })
}

mstep_pdz <-  function(y, posterior){
  temp_result_4 <- apply(posterior,2,function(z){       # posterior also has (docs*terms) rows but the n. of columns equals the n of topics. Each col has to be multiplied
    z*vect_corpus_PLSA_2                                # element-wise multiplication. To be summed
  })
  temp_result_5 <- cbind(y, temp_result_4)
  colnames(temp_result_5) <- c("doc", colnames(temp_result_4))
  temp_result_6 <- aggregate(. ~doc, data=as.data.frame(temp_result_5), sum)               # kind of sumif thread: https://stackoverflow.com/questions/21607464/what-is-the-equivalent-of-the-sumif-function-in-r
  pdznumerator <- temp_result_6[,-1]                                                               # remove first column
  pdzdenominator <- apply(pdznumerator, 2, sum)
  apply(pdznumerator,1, function(x){
    x * 1/pdzdenominator
  })
}

mstep_pz <-  function(posterior){
  temp_result_7 <- apply(posterior,2,function(w){       # posterior also has (docs*terms) rows but the n. of columns equals the n of topics. Each col has to be multiplied
    w*vect_corpus_PLSA_2                                  # element-wise multiplication. To be summed
  })
  pznumerator <- apply(temp_result_7,2,sum)
  pzdenominator <- sum_corpus                           # sum_corpus is a constant
  pznumerator/pzdenominator
}


#### calculate probs ####
max_iter <- 50                           # Hofmann 2001 p 184: apparently 20-50 iterations are typically sufficient
parameter.iter <- parameter.init
for(i_loop in 1:max_iter){
  # update posterior
  posterior.iter2 <- estep2(for_idx_plsa_swapped$doc, for_idx_plsa_swapped$semant_unit, parameter.iter)
  posterior_index_rownames <- paste0("d",for_idx_plsa_swapped$doc, "_", "t", for_idx_plsa_swapped$semant_unit)
  rownames(posterior.iter2) <- posterior_index_rownames

  # swap doc with term in expand grid so that we have doc = 1, terms =1,2,.... Right now it's the other way out
  posterior_reindex_rownames <-  paste0("d",for_idx_plsa$doc, "_", "t", for_idx_plsa$semant_unit)
  re_idx_posterior <- match(posterior_reindex_rownames,posterior_index_rownames)
  posterior_reindex <- posterior.iter2[re_idx_posterior,]

  # update parameters one by one
  pwz.iter <- mstep_pwz(for_idx_plsa$semant_unit, posterior_reindex)
  colnames(pwz.iter) <- termnames
  pwz.iter <- t(pwz.iter)

  pdz.iter <- mstep_pdz(for_idx_plsa_swapped$doc, posterior.iter2)
  colnames(pdz.iter) <- docnames
  pdz.iter <- t(pdz.iter)

  pz.iter <- mstep_pz(posterior.iter)

  # assemble
  parameter.iter <- list(pwz.iter, pdz.iter, pz.iter)
}
