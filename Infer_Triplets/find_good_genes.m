function good_genes = find_good_genes(tfdata, iunique, IA_BC_Amax, IA_BC_Bmax,IA_BC_Cmax,IABC_Amin,IABC_Bmin,IABC_Cmin,IABCsame)
    idxs = 1:max(iunique);
    [loggenemeans,~] = get_logg(tfdata,iunique);
    
    plikely_thresh = 0.6; %threshold for consideration
    ncells = max(iunique); ngenes = size(M,1);
    combinations = combnk(1:ncells,3); %possible combinations of 3
    p1_sla = 0.5; p2a_sla = 1/3*(1-p1_sla); p2b_sla = p2a_sla; p2c_sla = p2a_sla; %prior probabilities for slave genes
    psla = [p1_sla p2a_sla p2b_sla p2c_sla]; %vector of prior probabilities for slaves

    min_expressed = 0.2;
    min_frac_expressed = 0.2;
    good_genes = zeros(ngenes,1);
    for icomb = 1:length(combinations)
        iii = combinations(icomb,:);
        sz = [0 0 0];
        for k = 1:3
            sz(k) = sum(iunique==iii(k));
        end
        if min(sz)>5
            pgi_sla_1 = IABCsame(:,icomb); %p(gi | i slave, 1 distribution)
            pgi_sla_2Amax = 2*IA_BC_Amax(:,icomb); pgi_sla_2Amax(pgi_sla_2Amax<0) = 0; %(this is unlikely: should all be > 0
            pgi_sla_2Bmax = 2*IA_BC_Bmax(:,icomb); pgi_sla_2Bmax(pgi_sla_2Bmax<0) = 0; % but I did find one floating point error)
            pgi_sla_2Cmax = 2*IA_BC_Cmax(:,icomb); pgi_sla_2Cmax(pgi_sla_2Cmax<0) = 0;

            pgi_mas_Tmin = 3*[IABC_Amin(:,icomb) IABC_Bmin(:,icomb) IABC_Cmin(:,icomb)];
            pgi_mas_T = (1/2)*[pgi_mas_Tmin(:,2)+pgi_mas_Tmin(:,3) pgi_mas_Tmin(:,1)+pgi_mas_Tmin(:,3) pgi_mas_Tmin(:,2)+pgi_mas_Tmin(:,1)];
            pgi_sla = psla(1)*pgi_sla_1 + psla(2)*pgi_sla_2Amax + psla(3)*pgi_sla_2Bmax + psla(4)*pgi_sla_2Cmax; % p(g_i | slave)
            tot_expressed = [sum(tfdata(:,iunique==idxs(iii(1)))>0,2) sum(tfdata(:,iunique==idxs(iii(2)))>0,2) sum(tfdata(:,iunique==idxs(iii(3)))>0,2)];
            tot_cells = [sum(iunique==idxs(iii(1)))*ones(size(pgi_sla)) sum(iunique==idxs(iii(2)))*ones(size(pgi_sla)) sum(iunique==idxs(iii(3)))*ones(size(pgi_sla))];
            frac = tot_expressed./tot_cells;
            expressed = (tot_expressed > min_expressed).*(frac > min_frac_expressed);
            pgi_mas_T = pgi_mas_T.*expressed.*((sum(expressed,2)>1)*ones(1,3));

            pgi_mas_0 = (1/3)*sum(pgi_mas_T,2); %p(g_i | master)
            pgi_mas_T = [pgi_mas_T pgi_mas_0];  

            %find p(T|{g}) as a function of p(a_i=1)/p(a_i=0)
            pT_g = zeros(length(allodds),4); %pT_g as a function of prior odds
            for i = 1:length(allodds)
                prior_odds = allodds(i); %p(i master)/p(i slave)
                odds = 1+ prior_odds*pgi_mas_T./(pgi_sla*ones(1,4));
                odds(odds~=odds) = 1;
                odds(odds==Inf) = 1;
                %pT_g(i,:) = prod(odds).*[pA pB pC p0]/dot([pA pB pC p0],prod(odds));
                sum_log_odds = sum(log(odds));
                prododds = exp(sum_log_odds - max(sum_log_odds)); %divide by max(prod(odds)) to avoid Inf
                pT_g(i,:) = prododds.*topology_prior/dot(topology_prior,prododds);
            end
            [tlikely_ind, plikely_ind] = find_Tlikely(pT_g);
            if and(plikely_ind>plikely_thresh,tlikely_ind~=0)
                [pasym_gi_T, psym_gi_T, pmarker_gi_T] = find_gene_probabilities(pgi_mas_T, pgi_sla, pgi_sym, tlikely_ind, psla, odds0, loggenemeans(:,idxs(iii)), sym_thresh);
                pasym_gi_T(pasym_gi_T~=pasym_gi_T) = 0;
                psym_gi_T(psym_gi_T~=psym_gi_T) = 0;
                pmarker_gi_T(pmarker_gi_T~=pmarker_gi_T) = 0;
                genestatus = zeros(size(psym_gi_T,1));
                genestatus(pasym_gi_T>p_thresh_asym) = 1;
                genestatus(psym_gi_T>p_thresh_sym) = 1;
                genestatus(pmar_gi_T>p_thresh_sym) = 1;
                good_genes = (good_genes + genestatus)>0;
            end
        end
    end
end

    
    
