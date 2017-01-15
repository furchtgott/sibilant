%compute integrals and calculate topologies

%inputs:
%tfdata --> log-transformed data. rows = genes. columns = replicates (several per cell type). 
%idxs --> which cell types are we interested in out of all of the cell types in tfdata?
%labs --> cell type names
%iunique --> what is the cell type identity of each column of tfdata?
%iunique has same length as the # of columns in tfdata, and the max value
%of iunique is the number of cell types.
%P --> prior probability distribution P(mu, sigma). dimensions nbins x nbins
%Mu, Sigma --> axes for P. dimensions nbins x nbins

Params = get_parameters();
nbins = Params.nbins; 
mumin = Params.mumin; mumax = Params.mumax;
sigmin = Params.sigmin; sigmax = Params.sigmax;
[loggenemeans, loggenestds] = calc_log_mean_std(tfdata, iunique);

if Params.use_kde == 1
    [bandwidth,density,Mu,Sigma]=kde2d([loggenemeans(:) loggenestds(:)],nbins,[mumin sigmin],[mumax sigmax]); %kernel density estimation. 
    % Reference: Botev. Z.I., Grotowski J.F and Kroese D. P. (2010). Kernel density estimation via diffusion. Annals of Statistics. 
    % Volume 38, Number 5, Pages 2916--2957
    density(density(:)<0) = 0; %get rid of numerical errors (sum is ~10^-8)
    dm=Mu(1,2)-Mu(1,1);
    ds=Sigma(2,1)-Sigma(1,1);
    P = density/(dm*ds*sum(density(:)));
else
    dm = (mumax-mumin)/nbins;
    ds = (sigmax-sigmin)/nbins;
    [Mu,Sigma]=meshgrid(mumin:dm:mumax,sigmin:ds:sigmax);
    density = ones(size(Mu));
    P = density/(dm*ds*sum(density(:)));
end

%start parallel cluster
pc = parcluster('local');
pc.JobStorageLocation = strcat('/scratch/furchtg/',getenv('SLURM_JOB_ID'));
matlabpool(pc,16)

tic;
Integrals = calculate_integrals(tfdata, idxs, iunique, P, Mu, Sigma, Params.cutoff); 
toc;

clear pc


ncells = size(Integrals.IAB,1); ngenes = size(Integrals.IAB,3);
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



