## Combine data function

combineData <- function(data1, data2, name1 = NULL, name2 = NULL, delim = "_"){
  # read in data
  exp_pre <- data1
  exp_post <- data2
  if(!is.null(name1)){
    colnames(exp_pre) <- paste(colnames(exp_pre), name1, sep = delim)
  }
  if(!is.null(name2)){
    colnames(exp_post) <- paste(colnames(exp_post), name2, sep = delim)
  }
  
  # find all genes
  allgene <- unique(c(rownames(exp_pre), rownames(exp_post)))
  missG_pre <- setdiff(allgene, rownames(exp_pre))
  missG_post <- setdiff(allgene, rownames(exp_post))
  
  # create matrix with 0 for missing genes
  miss_pre <- matrix(0, ncol = dim(exp_pre)[2], nrow = length(missG_pre))
  rownames(miss_pre) <- missG_pre
  colnames(miss_pre) <- colnames(exp_pre)
  miss_post <- matrix(0, ncol = dim(exp_post)[2], nrow = length(missG_post))
  rownames(miss_post) <- missG_post
  colnames(miss_post) <- colnames(exp_post)
  
  # bind data
  new_pre <- rbind(exp_pre, miss_pre)
  new_post <- rbind(exp_post, miss_post)
  new_pre <- new_pre[order(rownames(new_pre)),]
  new_post <- new_post[order(rownames(new_post)),]
  print(mean(rownames(new_pre) == rownames(new_post)))
  data_new <- cbind(new_pre, new_post)
  return(data_new)
}

