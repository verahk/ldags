test_that("user-defined score gives same result as optimization routine", {
  
  data <- replicate(3, 0:1)
  nlev <- rep(2, 3)
  levels <- lapply(nlev-1, seq.int, from = 0)
  j <- 1
  pa <- 2:3 
  lookup <- rlang:::new_environment()
  
  
  for (method in c("tree", "ldag", "part")) {
    # expected outcome 
    counts <- compute_freq_table(data, nlev, j, pa)
    struct_exp <- optimize_partition(counts, levels[pa], ess = 1, 
                                     method = method, regular = T, verbose = T)
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
  
  data <- replicate(3, sample(0:1, 10, T))
  nlev <- rep(2, 3)
  levels <- lapply(nlev-1, seq.int, from = 0)
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
