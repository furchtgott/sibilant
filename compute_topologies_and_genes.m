%compute integrals and calculate topologies

%inputs:
%tfdata --> log-transformed data. rows = genes. columns = replicates (several per cell type). 
%idxs --> which cell types are we interested in out of all of the cell types in tfdata?
%labs --> cell type names
%iunique --> what is the cell type identity of each column of tfdata?
%iunique has same length as the # of columns in tfdata, and the max value
%of iunique is the number of cell types.

[loggenemeans, loggenestds] = calc_log_mean_std(tfdata, iunique);
Integrals = calculate_integrals(tfdata, idxs, iunique, Params)

ncells = length(idxs); ngenes = size(tfdata, 1);
combinations = combnk(1:ncells,3); %possible combinations of 3
selection = combinations; %triplets to consider (could choose a subset based on a preliminary distance metric)

ncandid = zeros(size(combinations,1),1); tlikely = zeros(size(ncandid)); plikely = zeros(size(ncandid));
makeplot = 0;
%find most likely topologies
for j=1:size(combinations,1)
    icomb = j;
    iii = combinations(icomb,:);
    trip_data = loggenemeans(:,idxs(iii));
    [tlikely_ind, plikely_ind, ncandidates] = process_one_triplet(icomb, Integrals, Params, trip_data, makeplot, labs(iii));
    tlikely(j) = tlikely_ind;
    plikely(j) = plikely_ind;
    ncandid(j) = ncandidates;
end

%find set of marker and transition genes
good_genes = zeros(ngenes,1);
for i=1:size(combinations,1)
    if and(plikely(i)>plikely_thresh,tlikely(i)~=0)
        icomb = i;
        iii = combinations(icomb,:);
        trip_data = loggenemeans(:,idxs(iii));
        [prob_asym_topol, prob_marker_topol] = calculate_gene_probabilities(icomb,Integrals, Params,trip_data,tlikely(icomb));
        genestatus = find_classes(prob_asym_topol, prob_marker_topol, tlikely(icomb), Params);
        good_genes = (good_genes + abs(genestatus))>0;
    end
end



