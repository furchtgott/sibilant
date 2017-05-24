%compute integrals and calculate topologies

%inputs:
%tfdata --> log-transformed data. rows = genes. columns = replicates (several per cell type). 
%iunique --> what is the cell type identity of each column of tfdata?
%iunique has same length as the # of columns in tfdata, and the max value
%of iunique is the number of cell types. Should take values between  1 and
%max(iunique)
%labs --> cell type names 

[loggenemeans, loggenestds] = calc_log_mean_std(tfdata, iunique);
ncells = length(unique(iunique)); ngenes = size(tfdata, 1);

Params = get_parameters();

% 1. Compute numeric integrals
Integrals = calculate_integrals(tfdata, iunique, Params);


% 2. Find most likely topologies for each triplet
combinations = combnk(1:ncells,3); %possible combinations of 3
ncandid = zeros(size(combinations,1),1); tlikely = zeros(size(ncandid)); plikely = zeros(size(ncandid));
makeplot = 0; %make plots?

for j=1:size(combinations,1)
    icomb = j;
    iii = combinations(icomb,:);
    trip_data = loggenemeans(:,iii);
    [tlikely_ind, plikely_ind, ncandidates] = process_one_triplet(icomb, Integrals, Params, trip_data, makeplot, labs);
    tlikely(j) = tlikely_ind;
    plikely(j) = plikely_ind;
    ncandid(j) = ncandidates;
end

% 3. Find union of marker and transition genes over all triplets
good_genes = zeros(ngenes,1);
for i=1:size(combinations,1)
    % only consider triplets with unique non-null topology
    if and(plikely(i)>Params.plikely_thresh,ncandid(i)==1) 
        icomb = i;
        iii = combinations(icomb,:);
        trip_data = loggenemeans(:,iii);
        [prob_asym_topol, prob_marker_topol] = calculate_gene_probabilities(icomb,Integrals, Params,trip_data,tlikely(icomb));
        genestatus = find_classes(prob_asym_topol, prob_marker_topol, tlikely(icomb), trip_data, Params);
        good_genes = (good_genes + abs(genestatus))>0;
    end
end

% 4. Recluster in subspace of good_genes in order to compute new iunique.
% Then go back to step 1.


