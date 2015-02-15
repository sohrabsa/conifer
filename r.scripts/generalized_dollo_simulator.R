  # Dollo model forward simulator
  # input: 
  # 1. a tree in ape's PhyloModel, \tau
  # 2. rate matrix (and state space) \Q
  # 3. initial distribution of insertion point \bar{\nu} - uniform on the branches with a point mass on the root \Omega
  # 4. initial distribution for the CTMC
  library(ape)
  library(VGAM)
  library(lattice)
  
  # TODO: add absorbing state
    # won't progress through the tree
  # TODO: sample before insertion
  # TODO: save ancestral states
  # TODO: generate a better tree 
  # TODO: add nucleotide statistics
  # TODO: fix the problem with the unrooted tree
  # TODO: assume there is no complete delete
  # TODO: incorporate b to the state space
  # TODO: the above process shouldn't be able to delete
  # TODO: make it compatible with are state space
  # TODO: add the emission model
  # TODO: add distribution test on the nucleotide on each leaf at each site chaning the root
  # TODO: add test for above and below processs
  
  # Inputs
  runAbove <- TRUE
  runBelow <- TRUE
  shallPrint <- TRUE
  
  
  
  ## CTMC
  # states mapping
  # states are represented as numbers 1 to length(states)
  # Will simulate a CTMC until it passes the length len
  # input: rate matrix Q, legnth (time) duration
  # output: final state (an integer number)
  
  
  # Simulate wait time (get the corresponding row to the current state i, simulate with exp(-Q(i, i)),   )
  # maybe collect states
  sampleCTMC <- function(init.state, duration, rateMatrix) {
    time.passed <- 0
    current.state <- init.state
    states <- seq(colnames(rateMatrix))
    pprint(current.state)
    # absorbing state
    if (rateMatrix[current.state, current.state] == 0) return(current.state)
        
    while (time.passed < duration) {
      # simulate wait time
      wait.time <- rexp(1, -rateMatrix[current.state, current.state])
      
      # choose a new state (shouldn't be the same one)
      print('***')
      print(rateMatrix[current.state, -c(current.state)])
      current.state <- sample(states[states != current.state], 1, prob = rateMatrix[current.state, -c(current.state)])
      pprint(current.state)
      # absorbing state
      if (rateMatrix[current.state, current.state] == 0) return(current.state)
      
      # update time passed
      time.passed <- time.passed + wait.time
    }
    
    print(paste('final current.state = ', current.state))
    current.state
  }
  
  #[1] "(0,0)" "(0,1)" "(0,2)" "(0,3)" "(1,0)" "(1,1)"
  #[7] "(1,2)" "(1,3)" "(2,0)" "(2,1)" "(2,2)" "(2,3)"
  #[13] "(3,0)" "(3,1)" "(3,2)" "(3,3)"
  
  # preorder traverse of the tree and collect states
  # input: tree, initial.state
  # output: a list, with tuples as elemetns, where each tuple is (treeNode, state)
  treeCTMC.edge.point <- function(tree, initial.state, edge.index, point, rateMatrix) {
    tree <- reorder(tree)
    tree$edge.length[edge.index] <- tree$edge.length[edge.index] - point
  
    final.index <- edge.index
    desc <- getDescendants(tree, tree$edge[edge.index, 2])
    if (length(desc) != 0) {
      last.leaf <- max(desc[desc %in% seq(tree$tip.label)])
      final.index <- which(tree$edge[,2] == last.leaf)
    } 
    
    # handle leaf edges
    if (tree$edge[edge.index, 2] %in% seq(length(tree$tip.label)))
      final.index <- edge.index
    # handle the root
    if (edge.index == 1 & point == 0) {
      final.index <- nrow(tree$edge)
    }
    
    result <- list()
    result[tree$edge[edge.index,1]] <- initial.state
    
    pprint(final.index)
    for (i in edge.index:final.index) {
      result[tree$edge[i,2]] <- sampleCTMC(result[[tree$edge[i,1]]], tree$edge.length[i], rateMatrix)
    }
    
    # Replace null with first absorbing state
    deleteStateIndex <- indexForAbsorbingState(rateMatrix)
    
    for (i in seq(length(result))) {
      if (is.null(result[[i]])) {
        result[[i]] <- deleteStateIndex
      }
    }
    pprint(result)
    result
  }
  
  
  # preorder traverse of the tree and collect states
  # Assumes it's not possible to completely remove
  # input: tree, initial.state
  # output: a list, with tuples as elemetns, where each tuple is (treeNode, state)
  # find the state for the event point at the child end of the edge.index
  treeCTMC.edge.point.exclude <- function(tree, initial.state, edge.index, point, rateMatrix) {
    tree <- reorder(tree)
    tree$edge.length[edge.index] <- tree$edge.length[edge.index] - point
    
    final.index <- edge.index
    desc <- getDescendants(tree, tree$edge[edge.index, 2])
    if (length(desc) != 0) {
      last.leaf <- max(desc[desc %in% seq(tree$tip.label)])
      final.index <- which(tree$edge[,2] == last.leaf)
    } 
    
    # handle leaf edges
    if (tree$edge[edge.index, 2] %in% seq(length(tree$tip.label)))
      final.index <- edge.index
    # handle the root
    if (edge.index == 1 & point == 0) {
      return(NULL)
    }
    
    result <- list()
    result[tree$edge[1,1]] <- initial.state
    
    pprint(edge.index)
    pprint(final.index)
    edges <- (1:nrow(tree$edge))[-c(edge.index:final.index)]
    pprint(edges)
  
    for (i in edges) {
      result[tree$edge[i,2]] <- sampleCTMC(result[[tree$edge[i,1]]], tree$edge.length[i], rateMatrix)
    }
        
    # simulate state for the event point
    result[tree$edge[edge.index, 2]] <- sampleCTMC(result[[tree$edge[i,1]]], point, rateMatrix)
    
    pprint(result)
    result
  }
  
  # choose the edge, then the point on the edge (is equal to the round-robin, one pass sampling)
  # order the edges in pre-order
  pickPointOnEdge <- function(tree, mu) {
    tree <- reorder(tree)
    probs <- c(tree$edge.length, 1/mu)
    
    # choose a branch
    br <- sample(seq(probs), 1, prob=probs)
    if (br == length(probs)) {
      return(c(1, 0))
    } else {
      # pick a point on this edge
      point <- runif(1, 0, tree$edge.length[br])
      return(c(br, point))
    }
  }
  
  driver.flat <- function(tree, nSites, initialDistribution, rateMatrix, mu) {
    
    state.labels <- list(colnames(rateMatrix[[1]]), colnames(rateMatrix[[2]]))
    states <- seq(length(state.labels))
    
    # Sample multiple sites
    sites <- list()
    for (i in seq(nSites)) {
      edgePoint <- pickPointOnEdge(tree, mu)
      pprint(edgePoint)
      # set the state at the root (normal state for cancer)
      temp1 <- NULL
      
      if (runAbove & ( !(edgePoint[1] == 1 & edgePoint[2] == 0) | !runBelow) ) {
        if (edgePoint[1] == 1 & edgePoint[2] == 0) {
          edgePoint <- pickPointOnEdge(tree, INF)
        }
        temp1 <- treeCTMC.edge.point.exclude(tree, 1, edgePoint[1], edgePoint[2], rateMatrix[[1]])
        rareEventState <- temp1[tree$edge[edgePoint[1], 2]]
      } else {
        rareEventState <- sample(states[[1]], 1, prob=pi)
      }
      
      # TODO: change the second component of the state space to reflect rare event
      pprint(temp1)
      
      temp2 <- NULL
      if(runBelow)
        temp2 <- treeCTMC.edge.point(tree, rareEventState, edgePoint[1], edgePoint[2], rateMatrix[[2]])
      
      if (is.null(temp1)) 
        temp1 <- temp2
      else if (is.null(temp2))
        temp2 <- temp1
      
      # merge the two
      sites[[i]] <- list()
      for(j in seq(length(temp1))) {
        sites[[i]][j] <- ifelse(is.null(temp1[[j]]), temp2[[j]], temp1[[j]])
      }
      
      pprint(sites[[i]])
    }
    
    # Extract taxa
    nTaxa <- length(tree$tip)
    taxa <- list()
    for (i in seq(nTaxa)) {
      taxon <- list()
      for (j in seq(nSites)) {
        taxon[[j]] <- sites[[j]][[i]]
      }
      taxa[[i]] <- taxon
    }
    
    # augment state space
    state.labels[length(states) + 1] <- '-'
    
    # Taxa to fasta
    fasta <- taxaToFasta(taxa, state.labels)
    
    print(fasta)
    
    # save the tree
    write.tree(tree, "~/Desktop/originalTree.nwk")
  }
  
  xx.taxa <- NULL
  xx.state.labels <- NULL
  xx.table <- NULL
  xx.emission <- NULL
  
  driver.ctmc <- function(tree, nSites, initialDistribution, rateMatrix, mu) {
    print(rateMatrix[[1]])
    
    state.labels <- list(colnames(rateMatrix[[1]]), colnames(rateMatrix[[2]]))
    states <- list(seq(length(state.labels[[1]])), seq(length(state.labels[[2]])))
                             
    # Sample multiple sites
    sites <- list()
    for (i in seq(nSites)) {
      edgePoint <- pickPointOnEdge(tree, mu)
      pprint(edgePoint)
      # set the state at the root (normal state for cancer)
      temp1 <- NULL
     
      if (runAbove & ( !(edgePoint[1] == 1 & edgePoint[2] == 0) | !runBelow) ) {
        if (edgePoint[1] == 1 & edgePoint[2] == 0) {
          edgePoint <- pickPointOnEdge(tree, INF)
        }
        root.state <- sample(states[[1]], 1, prob=initial.distribution[[1]])
        temp1 <- treeCTMC.edge.point.exclude(tree, root.state, edgePoint[1], edgePoint[2], rateMatrix[[1]])
        rareEventState <- temp1[tree$edge[edgePoint[1], 2]]
      } else {
        pprint(states[[1]])
        pprint(initial.distribution[[1]])
        rareEventState <- sample(states[[1]], 1, prob=initial.distribution[[1]])
        pprint(rareEventState)
      }
      
      # TODO: change the second component of the state space to reflect rare event
      pprint(rareEventState)
      tempState <- stringToCTMCPair(state.labels[[1]][rareEventState[[1]]])
      pprint(tempState)
      rareEventState <- grep(generateCNCTMCStateString(tempState[1]-1, 1), state.labels[[2]])
      pprint(rareEventState)
      pprint(temp1)
      
      temp2 <- NULL
      if(runBelow)
        temp2 <- treeCTMC.edge.point(tree, rareEventState, edgePoint[1], edgePoint[2], rateMatrix[[2]])
      
      if (is.null(temp1)) 
        temp1 <- temp2
      else if (is.null(temp2))
        temp2 <- temp1
      
      # merge the two
      sites[[i]] <- list()
      for(j in seq(length(temp1))) {
          sites[[i]][j] <- ifelse(is.null(temp1[[j]]), temp2[[j]] + nrow(rateMatrix[[1]]), temp1[[j]])
      }
      
      pprint(sites[[i]])
    }
    
    # Merge state spaces
    state.labels <- c(colnames(rateMatrix[[1]]), colnames(rateMatrix[[2]]))
    states <- seq(state.labels)

    xx.state.labels <<- state.labels
    
    # Extract taxa
    nTaxa <- length(tree$tip)
    taxa <- list()
    for (i in seq(nTaxa)) {
      taxon <- list()
      for (j in seq(nSites)) {
        taxon[[j]] <- sites[[j]][[i]]
      }
      taxa[[i]] <- taxon
    }
    
    xx.taxa <<- taxa
    
    # extract csv version 
    tableVersion <- CNfromTaxa(taxa, state.labels, '~/Desktop/dollo.cn.csv')
    print(tableVersion)
  
    xx.table <<- tableVersion
    
    # generate emission
    tableEmission <- tableVersion
    for (i in seq(nrow(tableVersion))) {
      tempEmission <- emissions(tableVersion$ref_counts[i], tableVersion$alt_counts[i])
      tableEmission$ref_counts[i] <- tempEmission[1]
      tableEmission$alt_counts[i] <- tempEmission[2]
    }
    write.csv(tableEmission, '~/Desktop/dollo.emission.csv')
    print(tableEmission)
    pplot(tableEmission)
    xx.emission <<- tableEmission
    # save the tree
    write.tree(tree, "~/Desktop/originalTree.nwk")
  }
  
  xx.tree <- NULL
  # get the parent of the chosen node
  # itterate till you get to the parrent again, stop
  example <- function() {
    nSites <- 48
    nTaxa <- 4
    set.seed(4)
    x <- rtree(nTaxa,T)
    x$tip.label <- paste0('t', 1:nTaxa)
    
    # SHOULD simulated above rare event?
    runAbove <<- TRUE
    runBelow <<- TRUE
    
    #x <- xx.tree
    plot(x)
    nodelabels()
    tiplabels();
    x <- reorder(x)
    xx.tree <<- x
    
    mu <- .01
    rateMatrix <- deleteRateMatrix(mu, size = length(states))
    print(rateMatrix)
    initial.distribution <- rep(1/length(states), length(states))
    # put zero for the delete state
    if (rateMatrix[length(states), length(states)] == 0) {
      initial.distribution[length(states)] <- 0
    }
    
    # (tree, nSites, initialDistribution, rateMatrix , states, state.labels, mu) {
    driver(x, nSites, initial.distribution, rateMatrix, states, state.labels, mu)
  }
  
  example.ctmc <- function() {
    nSites <- 48
    nTaxa <- 4
    #set.seed(4)
    maxCopyNumber <- 3
    x <- rtree(nTaxa,T)
    x$tip.label <- paste0('t', 1:nTaxa)
    
    # SHOULD simulated above rare event?
    runAbove <<- TRUE
    runBelow <<- TRUE

    # Plot and reorder the tree
    plot(x)
    nodelabels()
    tiplabels();
    x <- reorder(x)
    xx.tree <<- x
    
    rateMatrix <- list(generateCNRateMatrix.above(maxCopyNumber), generateCNRateMatrix.below(maxCopyNumber))
    mu <- .1
    
    # initial distribution
    initial.distribution <- list()
    for (r in rateMatrix) {
      initial.distribution[[length(initial.distribution) + 1]] <- degenerateInitDistribution(r, '(2,0)')
    }
    
    # (tree, nSites, initialDistribution, rateMatrix, mu) {
    driver.ctmc(x, nSites, initial.distribution, rateMatrix, mu)
  }
  
  ## Utils  
  CNfromTaxa <- function(taxa, state.labels, outPath) {
    # "sample_id","site_id","ref_counts","alt_counts","cluster_id"
    result <- NULL
    ref_counts <- c()
    alt_counts <- c()
    sample_id <- rep(seq(length(taxa)), each = length(taxa[[1]]))
    site_id <- rep(seq(length(taxa[[1]])), times = length(taxa))
    cluster_id <- rep(1, length(taxa) * length(taxa[[1]]))
    for (i in seq(length(taxa))) {
      taxon <- taxa[[i]]
      for (state in taxon) {
        tempState <- stringToCTMCPair(state.labels[state])
        ref_counts <- c(ref_counts, tempState[1])
        alt_counts <- c(alt_counts, tempState[2])
      }
    }
    
    result <- data.frame(sample_id, site_id, ref_counts, alt_counts, cluster_id)
    write.csv(result, outPath)
    result
  }
  
  fastaFromTaxa <- function(taxa, state.labels, outPath) {
    fasta <- c()
    for (i in seq(length(taxa))) {
      taxon <- taxa[[i]]
      fasta <- c(fasta, paste0('>t', i))
      fasta <- c(fasta, paste(state.labels[unlist(taxon)], collapse=''))
    }
    
    fileConn <- file(outPath)
    writeLines(fasta, fileConn)
    close(fileConn)
    fasta
  }
  
  # Returns a vector of absorbing states (with length 1 if there's only one absorbing state)
  indexForAbsorbingState <- function(rateMatrix) {
    which(diag(rateMatrix) == 0)
  }
  
  # Set to zero the probability of getting into any absorbing state
  removeAbsorbingState <- function(rateMatrix) {  
    for(i in seq(nrow(rateMatrix))) {
      if (rateMatrix[i, i] == 0) {
        for(j in seq(nrow(rateMatrix))) {
          rateMatrix[j, i] <- 0
        }
      }
    }
    
    for(i in seq(nrow(rateMatrix))) {
      rateMatrix[i, i] <- -sum(rateMatrix[i, -c(i)])
    }
    
    return(rateMatrix)
  }
  
  degenerateInitDistribution <- function(rateMatrix, degenStateStr) {
    states <- colnames(rateMatrix)
    as.numeric(states == degenStateStr)
  }
  
  deathlessInitialDistribution <- function(rateMatrix) {
    temp.dist <- rep(1/ncol(rateMatrix), ncol(rateMatrix))
    for (i in nrow(r)) {
      if (r[i, i] == 0)
        temp.dist[i] <- 0
    }
    temp.dist
  }
  
  simpleRateMatrix <- function() {
    result <- matrix(c(-3, 1, 1, 1,   1, -3, 1, 1,   1, 1, -3, 1,    1, 1, 1, -3   ), 4, 4)
    
    colnames(result) <- c("A", "T", "C", "G")
    rownames(result) <- c("A", "T", "C", "G")

    result
  }
  
  deleteRateMatrix <- function(deleteRate=1, size=5) {
   
    result <- matrix(c(0, 1, 1, 1, deleteRate,   
                       1, 0, 1, 1, deleteRate,   
                       1, 1, 0, 1, deleteRate, 
                       1, 1, 1, 0, deleteRate,     
                       0, 0, 0, 0, 0), size, size, byrow=TRUE)
    for (i in seq(size))
      result[i, i] <- -sum(result[i,])
    
    return(result)
  }
  
  generateCNCTMCStateString <- function(i, j) {
    paste0('(', i, ',', j, ')')
  }
  
  stringToCTMCPair <- function(stateStr) {
    #stateStr <- "(1 , 2)"
    result <- gsub('\\(', '', stateStr)
    result <- gsub('\\)', '', result)
    result <- gsub(" ", "", result)
    result <- strsplit(result, ",")
    as.numeric(result[[1]])
  }
  
  
  # when a mutation hasn't happened
  # Should not be able to go to a complete delete state
  generateCNRateMatrix.above <- function(maxCN = 4) {
    states <- c()
    rateMatrix <- matrix(0, (maxCN + 1), (maxCN + 1))
    
    # here b corresponds to state -1, therefore there's no b in this simulator
    for (i in 0:maxCN) { 
      states <- c(states, generateCNCTMCStateString(i, 0))  
    }
  
    for (i in seq(length(states))) {
      pair <- stringToCTMCPair(states[i])
      for (m in list(c(-1, 0), c(1, 0))) {
        if (   (pair[1] + m[1]) >= 0 & 
                 (pair[1] + m[1]) <= maxCN & 
                 (pair[2] + m[2]) >= 0 & 
                 (pair[2] + m[2]) <= maxCN) {
          tempState <- generateCNCTMCStateString(pair[1] + m[1], pair[2] + m[2])
          # prevent back mutation
          if ( (pair[2] != 0 | m[2] == 0)  &  (pair[1] != 0 | m[1] == 0))
            rateMatrix[i, grep(tempState, states)] <- 1
          
          # prevent complete deletion
            rateMatrix[2, 1] <- 0
        }
      }
    }
    
    # diagonal entries
    for (i in seq(states))
      rateMatrix[i, i] <- -sum(rateMatrix[i,])
    
    colnames(rateMatrix) <- states
    rownames(rateMatrix) <- states
    
    rateMatrix
  }
  
  # after a mutation (rare event occured)
  generateCNRateMatrix.below <- function(maxCN = 4) {
    states <- c()
    rateMatrix <- matrix(0, (maxCN + 1)^2, (maxCN + 1)^2, byrow = TRUE)

    # here b corresponds to state -1, therefore there's no b in this simulator
    for (i in 0:maxCN) { 
      for (j in 0:maxCN) { 
        states <- c(states, generateCNCTMCStateString(i, j))  
      }
    }
    
    for (i in seq(length(states))) {
      pair <- stringToCTMCPair(states[i])
      for (m in list(c(-1, 0), c(1, 0), c(0, -1), c(0, 1))) {
          if (   (pair[1] + m[1]) >= 0 & 
                 (pair[1] + m[1]) <= maxCN & 
                 (pair[2] + m[2]) >= 0 & 
                 (pair[2] + m[2]) <= maxCN) {
            tempState <- generateCNCTMCStateString(pair[1] + m[1], pair[2] + m[2])
            # prevent back mutation
            if ( (pair[2] != 0 | m[2] == 0)  &  (pair[1] != 0 | m[1] == 0))
              rateMatrix[i, grep(tempState, states)] <- 1
          }
      }
    }
    
    # diagonal entries
    for (i in seq(states))
      rateMatrix[i, i] <- -sum(rateMatrix[i,])
    
    colnames(rateMatrix) <- states
    rownames(rateMatrix) <- states
   
    rateMatrix
  }
  
  pprint <- function(a) {
    if (shallPrint) {
      arg <- deparse(substitute(a))
      print(paste0(arg, ' = ', a))
    }
  }
   
  getDescendants<-function(tree,node,curr=NULL){
    if(is.null(curr)) curr<-vector()
    daughters<-tree$edge[which(tree$edge[,1]==node),2]
    curr<-c(curr,daughters)
    w<-which(daughters>=length(tree$tip))
    if(length(w)>0) for(i in 1:length(w)) 
      curr<-getDescendants(tree,daughters[w[i]],curr)
    return(curr)
  }
  
  ### CTMC specific
  emission <- function(i, j) {
    poiMean <- 2000 
    trials <- rpois(1, poiMean)
    mean <- generateXi(i, j)
    s <- 1 # betaBinomialprecision
   
    mu <- mean
    rho <- 1/(1+s)
    major <- rbetabinom(1, trials, mu, rho)     
    
    minor <- trials - major
    c(major, minor)
  }
  
  muRhoFromMeanPrecision <- function(mean, s) {
    alpha <- mean * s
    beta <- (1-mean) * s
    mu <- alpha/(alpha + beta)
    rho <- 1/(1+ alpha+beta)
    return(c(mu, rho))
  }
  
  generateBetaBinomial <- function(n, mean, precision) {
    alpha <- mean * precision
    beta <- (1-mean) * precision
    p <- rbeta(1, alpha, beta)
    beta.binomial.sample <- rbinom(1, n, prob = p)
  }
  
  generateXi <- function(A, a) {
    DELTA <-  0.1
    xi <- 0
    if (A == 0 & a == 0)
      xi <- 0.5
    else if (A == 0)
      xi <- DELTA / (a + DELTA)
    else if (a == 0)
      xi <- A / (A + DELTA)
    else
      xi <- A / (A + a)
    return(xi)
  }
  
  pplot <- function(df) {
    # "sample_id","site_id","ref_counts","alt_counts","cluster_id"
    xyplot(ref_counts + alt_counts ~ site_id, df, auto.key = T)
  }
  
  #### Run examles 

  example.ctmc()
