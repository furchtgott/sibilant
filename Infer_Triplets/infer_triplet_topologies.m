cd('/n/home06/smelton1/Bayesian_Functions');

base_dir = '/n/home06/smelton1/Grun/';

load('/n/home06/smelton1/Grun/lgr_fulldata.mat');

cutoff = 0;
initial_cluster_number = 25;
nbins = 100; 
mumin = 0; mumax = 1.0;
sigmin = 0; sigmax = 0.75;
dm = -(mumin-mumax)/nbins;
ds = -(sigmin-sigmax)/nbins;
[Mu,Sigma]=meshgrid(mumin:dm:mumax,sigmin:ds:sigmax);
density = ones(size(Mu));
P = density/(dm*ds*sum(density(:)));

pc = parcluster('local');
pc.JobStorageLocation = strcat(base_dir,'Matlab_Storage');
parpool(pc,8)

clust_opts = statset('UseParallel',1);
par_kmed = @(X,K)(kmedoids(X,K,'options',clust_opts));
good_genes = ones(size(M,1),1);
dn_good_genes =[];

counter = 0;
go = counter==0;
while go
    spectral_M = spectral_mod(M(find(good_genes),:));
    iunique = par_kmed(spectral_M,initial_cluster_number);
    [~,korder] = sort(iunique);
    cell_cell = cov(M(find(good_genes),:));
    figure;
    imagesc(cell_cell(korder,korder));
    fn = strcat(base_dir,'cell_cell',int2str(counter),'.png');
    saveas(gcf,fn);
    close all;
    
    idxs = 1:max(iunique);
	tic;
	[~, ~, IABCsame, IABC_Amin, IABC_Bmin, IABC_Cmin, IA_BC_Amax, IA_BC_Bmax, IA_BC_Cmax] = calculate_integrals_both_parallel_test4(M, idxs, iunique, P, Mu, Sigma, cutoff); 
	toc;
    
    good_genes_new = find_good_genes(M, iunique, IA_BC_Amax, IA_BC_Bmax,IA_BC_Cmax,IABC_Amin,IABC_Bmin,IABC_Cmin,IABCsame);
    dn_good_genes = [dn_good_genes sum(good_genes(good_genes_new==0))];
    go = dn_good_genes(end)>140;
    good_genes = and(good_genes>0,good_genes_new>0);
    counter = counter + 1;
end

gg_spectral_M = spectral_mod(M(find(good_genes),:));
eva = evalclusters(gg_spectral_M,par_kmed,'gap','KList',[5:15],'B',50);

iunique = par_kmed(spectral_M,eva.OptimalK);
[~,korder] = sort(iunique);
cell_cell = cov(M(find(good_genes),:));
figure;
imagesc(cell_cell(korder,korder));
fn = strcat(base_dir,'cell_cell',int2str(counter),'.png');
saveas(gcf,fn);
close all;

idxs = 1:max(iunique);
tic;
[~, ~, IABCsame, IABC_Amin, IABC_Bmin, IABC_Cmin, IA_BC_Amax, IA_BC_Bmax, IA_BC_Cmax] = calculate_integrals_both_parallel_test4(M, idxs, iunique, P, Mu, Sigma, cutoff); 
toc;

ncells = max(iunique); ngenes = size(M,1);
combinations = combnk(1:ncells,3); %possible combinations of 3
p1_sla = 0.5; p2a_sla = 1/3*(1-p1_sla); p2b_sla = p2a_sla; p2c_sla = p2a_sla; %prior probabilities for slave genes
psla = [p1_sla p2a_sla p2b_sla p2c_sla]; %vector of prior probabilities for slaves

min_expressed = 0.2;
min_frac_expressed = 0.2;

[inferred_topologies, probs] = calc_top(combinations,psla, IABCsame,IA_BC_Amax,IA_BC_Bmax,IA_BC_Cmax, IABC_Amin, IABC_Bmin, IABC_Cmin, M, min_expressed, min_frac_expressed, idxs, iunique);

save(strcat(base_dir,'iter_oud_ints.mat'));
