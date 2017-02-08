function status = find_classes(pasym_gi_T, pmar_gi_T, tlik, data, Params)
    N_thresh_asym = Params.N_thresh_asym; 
    N_thresh_sym = Params.N_thresh_sym; 
    p_thresh_asym = Params.p_thresh_asym;
    p_thresh_sym = Params.p_thresh_sym;
    ngenes = size(pasym_gi_T,1);
    status = zeros(ngenes,1);
    % asym
    asymgenes = find(pasym_gi_T > p_thresh_asym); %find genes with prob(beta_i | g_i, T) > p_thresh_asym
    [~,asym_order] = sort(pasym_gi_T(asymgenes),'descend');
    asymgenes = asymgenes(asym_order(1:min(N_thresh_asym,length(asymgenes)))); %keep top N_thresh_asym of these genes
    [~,imax] = max(data(asymgenes,setdiff(1:3,tlik)),[],2); %sort into 2 classes depending on expression
    [~,group] = histc(imax,1:2);
    status(asymgenes(group==1)) = 1; %mark status in output vector: +/- 1
    status(asymgenes(group==2)) = -1;
    
    % markers for root cell type
    psym_gi_T = pmar_gi_T(:,tlik);
    symgenes = find(psym_gi_T > p_thresh_sym); %find genes with prob(alpha_i = 1| g_i, T) > p_thresh_sym
    [~,sym_order] = sort(psym_gi_T(symgenes),'descend');
    symgenes = symgenes(sym_order(1:min(N_thresh_sym,length(symgenes)))); %keep top N_thresh_sym of these genes
    status(symgenes) = 2; %mark status in output vector: +2
    
    % markers for leaf cell types
    daught = setdiff(1:3,tlik);
    daught1 = daught(1); % code: +3
    daught2 = daught(2); % code: -3
    markergenes1 = find(pmar_gi_T(:,daught1) > p_thresh_sym); %find genes with prob(alpha_i = 1| g_i, T) > p_thresh_sym
    markergenes2 = find(pmar_gi_T(:,daught2) > p_thresh_sym); %find genes with prob(alpha_i = 1| g_i, T) > p_thresh_sym
    [~,mar1_order] = sort(pmar_gi_T(markergenes1, daught1),'descend');
    markergenes1 = markergenes1(mar1_order(1:min(N_thresh_sym,length(markergenes1)))); %keep top N_thresh_sym of these genes
    [~,mar2_order] = sort(pmar_gi_T(markergenes2, daught2),'descend');
    markergenes2 = markergenes2(mar2_order(1:min(N_thresh_sym,length(markergenes2)))); %keep top N_thresh_sym of these genes
    status(markergenes1) = 3;
    status(markergenes2) = -3;
end