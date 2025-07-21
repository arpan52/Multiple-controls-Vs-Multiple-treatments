# Multiple-controls-Vs-Multiple-treatments

# Description
This repository contains all the original **R** codes used for the numerical studies in the paper entitled **On power optimal designs for multiple treatment comparisons** authored by **_Arpan Singh, Satya Prakash Singh_**.

# Instructions
 To reproduce the results, use the original **R** codes that are available in the **Main** branch. 

1) The following **R** packages are required to use the codes,\
      a) **nloptr**\
      b) **mvtnorm**\
      c) **tictoc**\
      d) **igraph**\
      e) **gtools**
      
3) Names of the main files and the **R** codes are self-explainatory. For example the file **Empirical power (Table 4)** contains all the **R** codes used to produce all the **Empirical powers** of all the designs given in the **Table 4**. 

4) Similarly, the **R** code **empirical_power_IUT_K6_1424343536.r**  produce the **empirical power** for **IUT** test associated to a bipartite graph with **K=6** comparison groups. Further, the set of pairs to be compared are **{(1,4), (2,4), (3,4), (3,5), (3,6)}**.

5) The **R** script **Form_of_MM_design_UIT_general_bipartite_graph_Final_1.R** implements **Theorem 3.5** from the article to determine the structure of the max–min optimal design for AN–power under a **general bipartite experimental graph**. The input to the code is a list of pairwise comparisons representing the edges of the experimental graph defined over K experimental groups. The output is a K-dimensional vector, where identical entries indicate the groups that should receive equal allocation under the max–min design. For example, the input edges1 <- list(c(1,2), c(1,3), c(1,4)) specifies a bipartite graph with edges (1,2), (1,3), and (1,4). The corresponding output is the vector (1, 2, 2, 2), indicating that the second, third, and fourth experimental groups should receive equal allocation, distinct from the first group.
