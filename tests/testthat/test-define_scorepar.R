test_that("user-defined score gives same result as optimization routine", {
  
  levels <- rep(list(0:1), 3)
  nlev <- lengths(levels)
  data <- sapply(levels, sample, size = 10, replace = T)
  j <- 1
  pa <- 2:3 
  lookup <- rlang:::new_environment()
  
  
  for (method in c("tree", "ldag", "part")) {
    # expected outcome 
    counts <- compute_freq_table(data, nlev, j, pa)
    struct_exp <- optimize_partition(counts, levels[pa], ess = 1, 
                                     method = method, regular = T, verbose = F)
    struct_exp <- struct_exp
    score_exp  <- sum(struct_exp$scores)
    
    # no lookup 
    scorepar <- define_scorepar(data, nlev, ess = 1, edgepf = 1, 
                                local_struct = method, regular = TRUE, lookup = NULL)
    score  <- BiDAG:::usrDAGcorescore(j, pa, length(nlev), scorepar)
    expect_equal(score, score_exp)
    
    # with lookup 
    id <- paste0(c(j, pa), collapse = "|")
    scorepar <- define_scorepar(data, nlev, ess = 1, edgepf = 1,
                                local_struct = method, regular = TRUE, lookup = lookup)
    score  <- BiDAG:::usrDAGcorescore(j, pa, length(nlev), scorepar)
    ## check score
    expect_equal(score, score_exp, ignore_attr = T)
    ## check score attribute
    expect_equal(attr(score, "struct"), struct_exp)     
    
    ## check score stored in lookup 
    expect_equal(lookup[[method]][[id]], score) 
  }
})


test_that("user-defined score gives correct score with none or a single parent", {
  
  levels <- rep(list(0:1), 3)
  nlev <- lengths(levels)
  data <- sapply(levels, sample, size = 10, replace = T)
  j <- 1
  lookup <- rlang:::new_environment()
  
  scorepar <- define_scorepar(data, nlev, ess = 1, edgepf = 1, 
                              local_struct = "tree", regular = TRUE, lookup = lookup)

  # no parents 
  pa <- integer(0)
  counts <- compute_freq_table(data, nlev, j, pa)
  score_exp  <- famscore_bdeu_1row(counts, 1, r = 2, q = 1, s = 1)
  score  <- BiDAG:::usrDAGcorescore(j, pa, length(nlev), scorepar)
  expect_equal(score, score_exp) 
  
  # one parent
  pa <- 2
  counts <- compute_freq_table(data, nlev, j, pa)
  score_exp  <- sum(famscore_bdeu_byrow(counts, 1, r = 2, q = 2, s = 1))
  score  <- BiDAG:::usrDAGcorescore(j, pa, length(nlev), scorepar)
  expect_equal(score, score_exp) 
  
})

test_that("user-defined score add correct edge penalty", {
  levels <- rep(list(0:1), 3)
  nlev <- lengths(levels)
  data <- sapply(levels, sample, size = 10, replace = T)
  
  edgepf <- 2
  scorepar <- define_scorepar(data, nlev, ess = 1, edgepf = edgepf, 
                              local_struct = "tree", regular = FALSE, lookup = NULL)
  
  j  <- 1
  pa <- 2
  res <- BiDAG:::usrDAGcorescore(j, pa, length(nlev), scorepar)
  
  # expected value
  famscore  <- compute_local_bdeu_score(data, levels, nlev, j, pa, ess = 1, method = "tree")
  res_exp   <- c(famscore) - log(edgepf)*length(pa)
  expect_equal(res, res_exp)
  
  j  <- 1
  pa <- 2:3
  res <- BiDAG:::usrDAGcorescore(j, pa, length(nlev), scorepar)
  
  # expected value
  famscore  <- compute_local_bdeu_score(data, levels, nlev, j, pa, ess = 1, method = "tree")
  res_exp   <- c(famscore) - log(edgepf)*length(pa)
  expect_equal(res, res_exp)
  
  
  
})
