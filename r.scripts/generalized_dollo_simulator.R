  # Dollo model forward simulator
  # input: 
  # 1. a tree in ape's PhyloModel, \tau
  # 2. rate matrix (and state space) \Q
  # 3. initial distribution of insertion point \bar{\nu} - uniform on the branches with a point mass on the root \Omega
  # 4. initial distribution for the CTMC
  library(ape)
  
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
  states <- NULL
  mu <- NULL
  Q <- NULL
  pi <- NULL
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
    pprint(current.state)
    # absorbing state
    if (rateMatrix[current.state, current.state] == 0) return(current.state)
        
    while (time.passed < duration) {
      # simulate wait time
      wait.time <- rexp(1, -rateMatrix[current.state, current.state])
      
      # choose a new state (shouldn't be the same one)
      current.state <- sample(states[states != current.state], 1, prob = rateMatrix[current.state, -c(current.state)])
      
      # absorbing state
      if (rateMatrix[current.state, current.state] == 0) return(current.state)
      
      # update time passed
      time.passed <- time.passed + wait.time
    }
    
    current.state
  }
  
  
  # preorder traverse of the tree and collect states
  # input: tree, initial.state
  # output: a list, with tuples as elemetns, where each tuple is (treeNode, state)
  treeCTMC.edge.point <- function(tree, initial.state, edge.index, point) {
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
      result[tree$edge[i,2]] <- sampleCTMC(result[[tree$edge[i,1]]], tree$edge.length[i], Q)
    }
    
    # change null with length(states) + 1
    for (i in seq(length(result))) {
      if (is.null(result[[i]])) {
        result[[i]] <- length(states) + 1
      }
    }
    pprint(result)
    result
  }
  
  
  # preorder traverse of the tree and collect states
  # It shouldne't be able to complete remove
  # input: tree, initial.state
  # output: a list, with tuples as elemetns, where each tuple is (treeNode, state)
  # find the state for the event point at the child end of the edge.index
  treeCTMC.edge.point.exclude <- function(tree, initial.state, edge.index, point) {
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
    
    # hack, temporaily change the rate matrix to remove the delte state
    qTemp <- removeAbsorbingState(Q)
  
    for (i in edges) {
      result[tree$edge[i,2]] <- sampleCTMC(result[[tree$edge[i,1]]], tree$edge.length[i], qTemp)
    }
        
    # simulate state for the event point
    result[tree$edge[edge.index, 2]] <- sampleCTMC(result[[tree$edge[i,1]]], point, qTemp)
    
    pprint(result)
    result
  }
  
  # remove the absorbing
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
  
  driver <- function(tree, nSites, initialDistribution, rateMatrix, states, state.labels, mu) {
    Q <<- rateMatrix
    pi <<- initialDistribution
    mu <<- mu
    states <<- states
    
    # Sample multiple sites
    sites <- list()
    for (i in seq(nSites)) {
      edgePoint <- pickPointOnEdge(tree, mu)
      pprint(edgePoint)
      # set the state at the root (normal state for cancer)
      temp1 <- NULL
      pprint(runAbove)
      pprint( !(edgePoint[1] == 1 & edgePoint[2] == 0 & !runBelow))
      
      if (runAbove & ( !(edgePoint[1] == 1 & edgePoint[2] == 0) | !runBelow) ) {
        if (edgePoint[1] == 1 & edgePoint[2] == 0) {
          edgePoint <- pickPointOnEdge(tree, INF)
        }
        temp1 <- treeCTMC.edge.point.exclude(tree, 1, edgePoint[1], edgePoint[2])
        rareEventState <- temp1[tree$edge[edgePoint[1], 2]]
      } else {
        rareEventState <- sample(states, 1, prob=pi)
      }
      
      pprint(temp1)
      
      temp2 <- NULL
      if(runBelow)
        temp2 <- treeCTMC.edge.point(tree, rareEventState, edgePoint[1], edgePoint[2])
      
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
    fasta <- c()
    for (i in seq(length(taxa))) {
      taxon <- taxa[[i]]
      fasta <- c(fasta, paste0('>t', i))
      fasta <- c(fasta, paste(state.labels[unlist(taxon)], collapse=''))
    }
    
    fileConn <- file("~/Desktop/dollo.fasta")
    writeLines(fasta, fileConn)
    close(fileConn)
    
    print(fasta)
    
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
    
    state.labels <- c("A", "T", "C", "G", '-')
    states <- seq(state.labels)
    mu <- .01
    #rateMatrix <- matrix(c(-3, 1, 1, 1,   1, -3, 1, 1,   1, 1, -3, 1,    1, 1, 1, -3   ), 4, 4)
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
    set.seed(4)
    maxCopyNumber <- 3
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
    
    rateMatrix.above <- generateCNRateMatrix.above(maxCopyNumber)
    rateMatrix.below <- generateCNRateMatrix.below(maxCopyNumber)
    
    state.labels <- colnames(rateMatrix.above)
    states <- seq(state.labels)
    mu <- .01
    print(rateMatrix.above)
    print(rateMatrix.below)
    
    initial.distribution <- rep(1/length(states), length(states))
    # put zero for the delete state
    if (rateMatrix[length(states), length(states)] == 0) {
      initial.distribution[length(states)] <- 0
    }
    
    # (tree, nSites, initialDistribution, rateMatrix , states, state.labels, mu) {
    driver(x, nSites, initial.distribution, rateMatrix, states, state.labels, mu)
  }
  
  ## Utils
  simpleRateMatrix <- function() {
    matrix(c(-3, 1, 1, 1,   1, -3, 1, 1,   1, 1, -3, 1,    1, 1, 1, -3   ), 4, 4)
  }
  
  deleteRateMatrix <- function(deleteRate=1, size=5) {
    pprint(size)
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
    
    print(states)
    
    for (i in seq(length(states))) {
      pair <- stringToCTMCPair(states[i])
      print(pair)
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
  
  #### Run examles 

  example.ctmc()
