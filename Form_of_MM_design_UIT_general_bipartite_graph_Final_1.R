rm(list = ls())

library(igraph)   # load for graphs
library(gtools)  # Load gtools for permutations


# Define the edges as pairs
#edges1 <- list(c(1,3), c(1,4), c(1,5), c(1,6), c(1,7), c(1,8), c(1,9), c(1,10), c(2,6), c(2,7), c(2,8), c(2,9), c(2,10))
#edges1 <- list(c(1,4), c(2,4), c(3,4), c(3,5), c(3,6))
#edges1 <- list(c(1,7), c(2,7), c(2,8), c(3,7), c(3,9), c(4,9), c(4,10), c(5,10), c(5,11), c(6,10))
#edges1 <- list(c(1,6), c(2,6), c(2,7), c(3,6), c(3,8), c(4,8),c(4,9),c(4,10), c(5,9))
#edges1 <- list(c(1,2), c(2,3), c(3,4))
edges1 <- list(c(1,2), c(1,3), c(1,4))
#edges1 <- list(c(1,3), c(1,4), c(1,5),c(2,3), c(2,4))

# Lexicographical Order for the comparisons
edges1 <- lapply(edges1, sort)
edges <- edges1[order(sapply(edges1, `[`, 1), sapply(edges1, `[`, 2))]


# Control set K_1
control_vertices<- unique(sapply(edges, function(x) x[1]))

# Treatment set K_2
treatment_vertices<- unique(sapply(edges, function(x) x[2]))

total_vertices<- c(control_vertices,treatment_vertices)

# find the graph from edges
# Convert list to edge matrix
edge_matrix <- do.call(rbind, edges)

# Create an undirected graph
g <- graph_from_edgelist(edge_matrix, directed = FALSE)

#######################################################################################################################
# For Path and Tree Graphs
# Extract unique vertices in order
vertices <- unique(unlist(edges))

# Print optimal design if valid
if (all(sapply(seq_len(length(edges) - 1), function(i) edges[[i]][2] == edges[[i + 1]][1]))) {
  print("Path graph")
  
  print("Optimal design:")
  print(c(1:ceiling(max(vertices)/2), rev(1:floor(max(vertices)/2))))
  return(FALSE)
}  

# Find the common elements across all edges (intersection of all vectors)
common_vertices <- Reduce(intersect, edges)

# Check if there is at least one common vertex in all edges
if (length(common_vertices) > 0) {
  result_vector <- rep(2, max(vertices))
  result_vector[common_vertices] <- 1
  print("Tree graph")
  print("Optimal design:")
  print(result_vector)
  return(FALSE)
}

# Check for Tree graph
if (all(diff(vertices) == 1)) {
  print("Path Graph")
  print("Optimal design:")
  print(c(1,rep(2,)))
  return(FALSE)
}

# Function to relabel every vertex sequentially while keeping the structure
relabel_edges_sequentially <- function(edges_list) {
  # Extract all unique vertices in the order they appear
  unique_vertices <- unique(unlist(edges_list))
  
  # Create a mapping from old vertex labels to new sequential numbers
  relabel_map <- setNames(seq_along(unique_vertices), as.character(unique_vertices))
  
  # Apply relabeling to edges using the mapping, ensuring valid numeric output
  relabeled_edges <- lapply(edges_list, function(edge) c(relabel_map[as.character(edge[1])], relabel_map[as.character(edge[2])]))
  # Convert back to numeric
  relabeled_edges <- lapply(relabeled_edges, as.numeric)  
  # Convert the list to a matrix
  relabeled_edges_matrix <- do.call(rbind, relabeled_edges)
  
  return(relabeled_edges_matrix)
}


######################################################################################################################

# Function to check if two graphs are isomorphic for the case of odd number of edges
are_isomorphic_odd <- function(edges1, edges2,control_vertices_1,treatment_vertices_1,control_vertices_2,treatment_vertices_2){
  
  
  # Convert edge lists to relabeled matrices
  edge_matrix1 <- relabel_edges_sequentially(edges1)
  edge_matrix2 <- relabel_edges_sequentially(edges2)
  
  # Create graph objects
  g1 <- graph_from_edgelist(edge_matrix1, directed = FALSE)
  g2 <- graph_from_edgelist(edge_matrix2, directed = FALSE)
  
  # Check if both graphs are connected before proceeding
  if (ecount(g1) < 2 || ecount(g2) < 2 || !is_connected(g1) || !is_connected(g2)) {
    return(FALSE)
  }
  
  # Check for isomorphism using igraph
  if (igraph::isomorphic(g1, g2)) {
    
    
    for (m in 1:length(control_vertices_1)){
      
      if(degree(g)[total_vertices[control_vertices_1[m]]]!= degree(g)[treatment_vertices_2[length(control_vertices_1)-m+1]])
      { return(FALSE)  # Degree mismatch → Not isomorphic
      }}
    for (m in 1:length(control_vertices_2)){
      if(degree(g)[total_vertices[control_vertices_2[m]]]!= degree(g)[treatment_vertices_1[m]])
      { return(FALSE)  # Degree mismatch → Not isomorphic
      }}
    
    return(TRUE)
  }
  
  return(FALSE)
}


# Code to execute if the length is odd

if (length(edges) %% 2 == 1) {
  edge_set1 <- edges[1:((length(edges) - 1) / 2)]
  edge_set2 <- edges[(((length(edges) - 1) / 2) + 2):length(edges)]
  
  
  control_vertices_1<- unique(sapply(edge_set1, function(x) x[1]))
  treatment_vertices_1<- unique(sapply(edge_set1, function(x) x[2]))
  control_vertices_2<- unique(sapply(edge_set2, function(x) x[1]))
  treatment_vertices_2<- unique(sapply(edge_set2, function(x) x[2]))
  
  if(are_isomorphic_odd(edge_set1, edge_set2,control_vertices_1,treatment_vertices_1,control_vertices_2,treatment_vertices_2)){
    
    for (m in 1:length(control_vertices_1)){
      total_vertices[control_vertices_1[m]]<- treatment_vertices_2[length(control_vertices_1)-m+1]
    }
    for (m in 1:length(control_vertices_2)){
      total_vertices[control_vertices_2[m]]<- treatment_vertices_1[m]
    }
    
    print("Optimal design:")
    print(total_vertices)
    return(FALSE)}
  
  
}

#########################################################################################################################

# Function to check if two graphs are isomorphic for the case of odd number of edges
are_isomorphic_even <- function(edges1, edges2,control_vertices_1,treatment_vertices_1,control_vertices_2,treatment_vertices_2){
  
  
  # Convert edge lists to relabeled matrices
  edge_matrix1 <- relabel_edges_sequentially(edges1)
  edge_matrix2 <- relabel_edges_sequentially(edges2)
  
  # Create graph objects
  g1 <- graph_from_edgelist(edge_matrix1, directed = FALSE)
  g2 <- graph_from_edgelist(edge_matrix2, directed = FALSE)
  
  # Check if both graphs are connected before proceeding
  if (ecount(g1) < 2 || ecount(g2) < 2 || !is_connected(g1) || !is_connected(g2)) {
    return(FALSE)
  }
  
  # Check for isomorphism using igraph
  if (igraph::isomorphic(g1, g2)) {
    
    for (m in 1:length(control_vertices_1)){
      
      if(degree(g)[total_vertices[control_vertices_1[m]]]!= degree(g)[control_vertices_2[length(control_vertices_1)-m+1]])
      { return(FALSE)  # Degree mismatch → Not isomorphic
      }}
    for (m in 1:(length(treatment_vertices_1)-1)){
      if(degree(g)[ total_vertices[treatment_vertices_1[m]]]!= degree(g)[treatment_vertices_2[m+1]])
      { return(FALSE)  # Degree mismatch → Not isomorphic
      }}
        return(TRUE)
  }
  
  return(FALSE)
}

# Code to execute if the length is even

if (length(edges) %% 2 == 0) {
  edge_set1 <- edges[1:((length(edges)) / 2)]
  edge_set2 <- edges[(((length(edges)) / 2) + 1):length(edges)]
  
  
  control_vertices_1<- unique(sapply(edge_set1, function(x) x[1]))
  treatment_vertices_1<- unique(sapply(edge_set1, function(x) x[2]))
  control_vertices_2<- unique(sapply(edge_set2, function(x) x[1]))
  treatment_vertices_2<- unique(sapply(edge_set2, function(x) x[2]))
  
  if(are_isomorphic_even(edge_set1, edge_set2,control_vertices_1,treatment_vertices_1,control_vertices_2,treatment_vertices_2)){
    
    
    for (m in 1:length(control_vertices_1)){
      total_vertices[control_vertices_1[m]]<- control_vertices_2[length(control_vertices_1)-m+1]
    }
    for (m in 1:(length(treatment_vertices_1)-1)){
      total_vertices[treatment_vertices_1[m]]<- treatment_vertices_2[m+1]
    }
    
    print("Optimal design:")
    print(total_vertices)
    return(FALSE)
  }
  
}


#######################################################################################################################
# For a general Bipartite graph which can not be written in two isomorphic subgraphs


sets <- list(control_vertices, treatment_vertices)

for (s in sets) {
  for (i in s) {
    for (j in Filter(function(x) degree(g)[x] == degree(g)[i], s[s > i])) {  # Select elements greater than i in s
      
      print(c(i, j))
      
      
      # Find vertices connected to i
      V_i  <- unique(unlist(lapply(edges, function(x) {
        if (x[1] == i) return(x[2])
        else if (x[2] == i) return(x[1])
        else return(NULL)  # This avoids adding NULL
      })))
      
      # Find vertices connected to j
      V_j <- unique(unlist(lapply(edges, function(x) {
        if (x[1] == j) return(x[2])
        else if (x[2] == j) return(x[1])
        else return(NULL)  # This avoids adding NULL
      })))
      

      
      # check the first part of the theorem    
      if (setequal(V_i, V_j)) {total_vertices[total_vertices == i] <- j
      print('Optimal Design')
      print(total_vertices)  
      next}
      
      # Find vertices connected to i but not to j
      U_ij <- setdiff(V_i, V_j)
      # Find vertices connected to j but not to i
      U_ji <- setdiff(V_j, V_i)
      
      # Find vertices connected to both i and j
      W_ij <- intersect(V_i, V_j)
      # Find vertices connected to i or j
      U<-  union(U_ij,U_ji)
      
      ##########################################################################################################
      if (all(c(i, j) %in% treatment_vertices))
      {
        # find set I_ij and J_ij
        I_K1_ij<-c()
        I_K2_ij<-c()
        
        
        # Define the set of nodes to check the intersection
        check_set <- union(c(i,j),W_ij)
        
        # find set I_K1_ij
        for (k in setdiff(control_vertices, W_ij)){
          for (k_u in U){
            # find all paths from k to k_u
            paths <- all_simple_paths(g, from = k, to = k_u)
            
            # If at least one path does NOT have anything in common with check_set
            if (any(!sapply(paths, function(path) any(path %in% check_set)))) {
              I_K1_ij <- c(I_K1_ij, k)  # Add k to I_K1_ij
            }
            
          }
          
        }
        I_K1_ij <- unique(c(I_K1_ij, U))
        # find set I_K2_ij
        for (k in setdiff(treatment_vertices, c(i,j))){
          for (k_u in U){
            # find all paths from k to k_u
            paths <- all_simple_paths(g, from = k, to = k_u)
            
            # If at least one path does NOT have anything in common with check_set
            if (any(!sapply(paths, function(path) any(path %in% check_set)))) {
              I_K2_ij <- c(I_K2_ij, k)  # Add k to I_K2_ij
            }
            
          }
          
        }
        I_K2_ij <- unique(I_K2_ij)}
      ###################################################################################################################
      
      if (all(c(i, j) %in% control_vertices))
      {
        # find set I_ij and J_ij
        J_K1_ij<-c()
        J_K2_ij<-c()
        
        # Define the set of nodes to check the intersection
        check_set <- union(c(i,j),W_ij)
        
        # find set J_K1_ij
        for (k in setdiff(control_vertices, W_ij)){
          for (k_u in U){
            # find all paths from k to k_u
            paths <- all_simple_paths(g, from = k, to = k_u)
            
            # If at least one path does NOT have anything in common with check_set
            if (any(!sapply(paths, function(path) any(path %in% check_set)))) {
              J_K1_ij <- c(J_K1_ij, k)  # Add k to J_K1_ij
            }
          }
          
        }
        J_K1_ij <- unique(J_K1_ij)
        # find set J_K2_ij
        for (k in setdiff(treatment_vertices, c(i,j))){
          for (k_u in U){
            # find all paths from k to k_u
            paths <- all_simple_paths(g, from = k, to = k_u)
            
            # If at least one path does NOT have anything in common with check_set
            if (any(!sapply(paths, function(path) any(path %in% check_set)))) {
              J_K2_ij <- c(J_K2_ij, k)  # Add k to J_K2_ij
            }
            
          }
          
        }
        J_K2_ij <- unique(c(J_K2_ij, U))
      }
      ################################################################################################################
      # for i,j in K_2 set, that is the set of treatments
      if(all(c(i,j) %in% treatment_vertices)){
        if (is.null(I_K1_ij) || is.null(I_K2_ij) || length(I_K1_ij) == 0 || length(I_K2_ij) == 0) {
          print("No isomorphic subgraph possile")
          next  # Skip this iteration and move to the next i, j
        }
        matching_edges <- Filter(function(edge) {(edge[1] %in% I_K1_ij && edge[2] %in% I_K2_ij)}, edges)}
      
      # for i,j in K_1 set, that is the set of controls
      if(all(c(i,j) %in% control_vertices)){
        if (is.null(J_K1_ij) || is.null(J_K2_ij) || length(J_K1_ij) == 0 || length(J_K2_ij) == 0) {
          print("No isomorphic subgraph possile")
          next  # Skip this iteration and move to the next i, j
        }
        matching_edges <- Filter(function(edge) {(edge[1] %in% J_K1_ij && edge[2] %in% J_K2_ij)}, edges)}
      
      
      # Get all subsets of matching_edges
      subgraph_edges <- unlist(lapply(0:length(matching_edges), function(k) combn(matching_edges, k, simplify = FALSE)), recursive = FALSE)[-1]
      
      
      # Function to check if two subgraphs are disconnected
      are_disconnected <- function(sub1, sub2) {
        # Convert edges to character format for comparison
        sub1_edges <- apply(do.call(rbind, sub1), 1, paste, collapse = "-")
        sub2_edges <- apply(do.call(rbind, sub2), 1, paste, collapse = "-")
        
        # Check if there is any common edge
        disconnected <- length(intersect(sub1_edges, sub2_edges)) == 0
        return(disconnected)
      }
      
      
      # Function to check if two graphs are isomorphic and display edges if they are
      are_graphs_isomorphic <- function(edges1, edges2) {
        
        if (!are_disconnected(edges1, edges2)) {
          return(FALSE)
        }
        
        # Convert edge lists to relabeled matrices
        edge_matrix1 <- relabel_edges_sequentially(edges1)
        edge_matrix2 <- relabel_edges_sequentially(edges2)
        
        # Create graph objects
        g1 <- graph_from_edgelist(edge_matrix1, directed = FALSE)
        g2 <- graph_from_edgelist(edge_matrix2, directed = FALSE)
        
        # Check if both graphs are connected before proceeding
        if (ecount(g1) < 2 || ecount(g2) < 2 || !is_connected(g1) || !is_connected(g2)) {
          return(FALSE)
        }
        
        # Check for isomorphism using igraph
        if (igraph::isomorphic(g1, g2)) {
          
          # Check if corresponding vertices have the same degree in external graph g
          for (l in seq_along(edges1)) {
            
            if (degree(g)[edges1[[l]][1]] != degree(g)[edges2[[l]][1]] ||
                degree(g)[edges1[[l]][2]] != degree(g)[edges2[[l]][2]]) {
              
              return(FALSE)  # Degree mismatch → Not isomorphic
            }
          }
          
          return(TRUE)
        }
        
        return(FALSE)
      }
      
      # Find isomorphic pairs
      for (u in 1:(length(subgraph_edges) - 1)) {
        for (v in (u + 1):length(subgraph_edges)) {
          
          if(are_graphs_isomorphic(subgraph_edges[[u]], subgraph_edges[[v]]))
            
            
          {          # Flatten the edge lists
            e1 <- unlist(subgraph_edges[[u]])  # Extract nodes from Graph 1
            e2 <- unlist(subgraph_edges[[v]])  # Extract corresponding nodes from Graph 2
            
            
            
            # Replace i with j in total_vertices
            total_vertices[total_vertices == i] <- j
            
            # Update total_vertices based on mappings in e1 and e2
            for (m in seq_along(e1)) {
              total_vertices[total_vertices == e1[m]] <- e2[m]
            }
            
          }
        }
      }
      
    }
  }
}
print("Optimal design:")
print(total_vertices)