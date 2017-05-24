# sibilant

MATLAB source code for Furchtgott et al. paper, "Discovering sparse transcription factor codes for cell states and state transitions during development." 

(README and documentation under development)

The code consists of three main functions. These functions are called in the script compute_topologies_and_genes.m which presents a typical workflow.

1. calculate_integrals.m computes the numeric integrals based on cluster identities and gene expression data and returns a MATLAB structure with fields for each numeric integral.

2. process_one_triplet.m computes probabilities of topology for a triplet of cell types using the numeric integrals. 

3. compute_gene_probabilities.m computes probabilities of gene classes (marker, transition) for a triplet of cell types assuming a particular topology and using the numeric integrals. 


All parameters are set in the function get_parameters() which returns a MATLAB structure with fields for each parameter. Currently the actual code for this function must be edited in order to change the values of parameters; it would be much better for the future to have the function read in a text file with the parameter values. 

Given new data, where to start, which parameters to set, etc.?

A. input data:

tfdata --> log2-transformed gene expression data. rows = genes. columns = replicates (several per cell type). 

iunique --> vector providing the cell type/cluster identity of each column of tfdata. iunique has same length as the # of columns in tfdata, and the max value of iunique is the number of cell types (or clusters). Should take values between 1 and max(iunique). 

labs --> cell type (cluster) names 

For single-cell data, you need to start with a seed cluster identity for iunique. We have used spectral k-means clustering and spectral t-SNE (using the Seurat R package).

B. Running the code: parameters for calculate_integrals.m

The numeric integral calculation can be quite time-consuming, particularly for a large number of cell types (and the corresponding factorial number of triplets). There is an option (use_parallel and parallel_pool_size) for running the code in parallel using a parfor loop in a parallel computing environment. 

mumin, mumax, sigmin, sigmax define the bounds of the prior P(mu, sigma) of means and standard deviations of the distribution of gene expression in a cluster. These will need to be adjusted depending on the dataset. 

use_kde if set to 1 will use kde-estimated empirical prior; if set to 0 will use uniform (flat) prior

nbins is the number of bins for the numeric integration

cutoff is an expression threshold for numeric integration for the downregulation pattern. This sets a bound on the limits of integration. In the integrals where muA < muB and muA < muC, then this sets an additional limit muB > cutoff and muC > cutoff; i.e. max(cutoff, muA ) < muB and max(cutoff, muA) < muC. See code for details. 
cutoff should be set either to a reasonable threshold value for expression or to 0 in which case this feature will be ignored.


C. running the code: parameters for process_one_triplet.m

Once the numeric integrals are calculated, the next function computes probabilities of topologies given data for each triplet. The main parameter to modify here is “oddsrange”, i.e. the range to vary prior odds p(beta_i = 1)/p(beta_i = 0). See Figure 2—figure supplement 1 B and D in the paper. You want to be scanning a range in which (for most triplets): for low values of prior odds, p(T|{g}) is p(T) = 1/4 (no data), and for high values of prior odds, p(T|{g}) goes to the null hypothesis. 
1e-6 to 1e2 is a good starting point, but I would recommend running process_one_triplet for a few triplets with makeplot = 1 to see the plots and check the range.

The result of this function is: tlikely (most likely non-null topology over oddsrange), plikely (probability of this topology), ncandid (how many candidate non-null topologies above plikely_thresh)

D. running the code: parameters for compute_gene_probabilities.m and find_classes.m

For a given triplet we can use compute_gene_probabilities.m to calculate the probabilities of each gene being a transition gene (prob_asym_topol) or marker gene (prob_marker_topol) given a) data, b) inferred topology, c) a given value odds0 of prior odds of being transition gene p(beta_i = 1)/p(beta_i = 0). 

Absolute probabilities will depend on odds0 but not the ordering of genes. The value of odds0 should be chosen in light of oddsrange.

There is an option (filter_min = 1) to filter transition and marker genes by setting expression thresholds log_asym_min and log_sym_min respectively. This ensures that 
transition genes have expression above the expression threshold in at least 2 cell types and marker genes have expression above expression threshold for the “marked" cell type. If using this option then log_asym_min and log_sym_min should be set to reasonable values. 

Function find_classes.m returns a vector status with length ngenes. The entry for each gene reflects whether it is a transition or marker gene for that particular triplet: 
status = 1 or -1: transition gene downregulated in one or the other leaf cell types
status = 2: marker gene for root cell type
status = 3 or -3: marker gene for one or the other leaf cell types

To define classes, need to use a threshold for transition or marker gene probabilities (p_thresh_asym and p_thresh_sym respectively). We can also chose the top N_thresh_asym genes for the transition gene class or the the top N_thresh_sym genes for each of the three marker gene classes. 

Note that the number of transition/marker genes is going to depend on user-defined odds0 and probability thresholds. 

E. Running code for single-cell data

The union of the genes in different classes for different triplets can be used to recluster the cells (using t-SNE or spectral k-means etc.). After reclustering and obtaining a new set of cluster identities iunique(2?), it’s necessary to recompute the numeric integrals. 