function Params = get_parameters()

Params.use_kde = 0;
Params.nbins = 2^8; 
Params.mumin = 2; 
Params.mumax = 14;
Params.sigmin = 0; 
Params.sigmax = 0.75;
Params.cutoff = log2(100);
Params.psla = [1/2 1/6 1/6 1/6]; %vector of prior probabilities for non-transition genes
Params.topology_prior = [1/4 1/4 1/4 1/4];
Params.plikely_thresh = 0.5;
Params.filter_min = 0;%ignore genes with values under threshold?
Params.log_asym_min = log2(100); %threshold for expression
Params.log_sym_min = log2(50); %threshold for expression
Params.allodds = logspace(-4,2,50); %vary prior odds p(i master)/p(i slave)

%parameters for gene probabilities
Params.N_thresh_sym = 1000; 
Params.N_thresh_asym = 1000; %thresholds for defining gene classes
Params.p_thresh_asym = 0.8; 
Params.p_thresh_sym = 0.8; 
Params.odds0 = 0.05;

end