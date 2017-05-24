function Params = get_parameters()

%numeric integration parameters
Params.use_kde = 0; %1 -> use kde-estimated empirical prior;  0 -> use uniform prior
Params.use_parallel = 0; %compute integrals in parallel?
Params.parallel_pool_size = 16; %number of processors for parallelization
Params.nbins = 2^8; %number of bins for numeric integration
Params.mumin = 2; % min and max for mean prior
Params.mumax = 14;
Params.sigmin = 0; % min and max for std prior
Params.sigmax = 0.75;
Params.cutoff = log2(100); % expression threshold for numeric integration

%parameters for inferring topologies
Params.psla = [1/2 1/6 1/6 1/6]; %vector of prior probabilities for non-transition genes
Params.topology_prior = [1/4 1/4 1/4 1/4]; %topology prior [p(A) p(B) p(C) p(0)]
Params.plikely_thresh = 0.6; % threshold for probability of topology given data
Params.oddsrange = logspace(-4,2,50); %range to vary prior odds p(beta_i = 1)/p(beta_i = 0)

%parameters for gene probabilities
Params.odds0 = 0.05; %transition gene prior for calculating gene probabilities
Params.filter_min = 0; %ignore genes with values under threshold?
Params.log_asym_min = log2(100); %threshold for transition gene expression
Params.log_sym_min = log2(50); %threshold for marker gene expression

%parameters for gene classes
Params.N_thresh_sym = 1000; % max number of genes in gene classes
Params.N_thresh_asym = 1000; 
Params.p_thresh_asym = 0.8; % probability thresholds for defining gene classes
Params.p_thresh_sym = 0.8; 

end