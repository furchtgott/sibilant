function [inferred_topologies, probs] = calc_top(combinations,psla, IABCsame,IA_BC_Amax,IA_BC_Bmax,IA_BC_Cmax, IABC_Amin, IABC_Bmin, IABC_Cmin, tfdata, min_expressed, min_frac_expressed, idxs, iunique)
    inferred_topologies = zeros(length(combinations),1);
    probs = zeros(length(combinations),1);
    topology_prior = [0.25 0.25 0.25 0.25]; % uniform prior over possible topologies (including null)
    allodds = logspace(-50,0,100);
    for icomb = 1:length(combinations)
        iii = combinations(icomb,:);
        
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
        probs(icomb) = plikely_ind;
        inferred_topologies(icomb) = tlikely_ind;
    end
end
