function [prob_asym_topol, prob_marker_topol] = calculate_gene_probabilities(icomb,Integrals, Params,trip_data,tlik)
    %See Equations (27) and (28)
    odds0 = Params.odds0;
    psla = Params.psla;
    [pgi_mas_T, pgi_sla, pgi_sym] = calc_pgi_topo(icomb,Integrals, Params, trip_data); %get probability of data given different hypotheses
    % transition genes
    odds_asym_topol = odds0*pgi_mas_T(:,tlik)./pgi_sla; % p(beta_i = 1 | g_i, T)/p(beta_i = 0 | g_i, T)
    prob_asym_topol = odds_asym_topol./(1+odds_asym_topol); % p(beta_i = 1 | g_i, T)
    prob_asym_topol(prob_asym_topol~=prob_asym_topol) = 0;  %get rid of NaNs  
 
    p_mas = odds0/(1+odds0); % prior p(beta_i=1)
    p_sla = 1/(1+odds0); % prior p(beta_i=0)
    pg_T = p_mas*pgi_mas_T(:,tlik) + p_sla*pgi_sla; % p(g_i | T)
    
    %markers
    %pg_T as above
    prob_marker_topol = zeros(size(pgi_sym));
    for t=1:3
        psym = pgi_sym(:,t); % p(g_i | T, gene i marker for cell type t) i.e. p(g_i | T, beta_i = 0, alpha_i = 1)
        pnotsym = pg_T-p_sla*psla(t+1)*psym; % p(g_i | T, gene i NOT marker for cell type t). psla(tlik+1) is p(mu_A > mu_BC | b_i = 0) (or equiv) 
        prior_sym_odds = (p_sla*psla(t+1))/(1-p_sla*psla(t+1)); % p(alpha_i = 1 for cell type t)/p(alpha_i = 0 for cell type t);
        odds_marker_topol = prior_sym_odds*psym./pnotsym; % p(alpha_i = 1 | g_i, T)/p(alpha_i = 0 | g_i, T)
        prob_marker_topol(:,t) = odds_marker_topol./(odds_marker_topol+1); % p(alpha_i = 1 for cell type t | g_i, T)
        prob_marker_topol(odds_marker_topol == Inf,t) = 1;
        prob_marker_topol(prob_marker_topol(:,t)~=prob_marker_topol(:,t),t) = 0;
    end

end

function [pgi_mas_T, pgi_sla, pgi_sym] = calc_pgi_topo(icomb,Integrals, Params, data)

    log_asym_min = Params.log_asym_min;
    log_sym_min = Params.log_sym_min;
    filter_min = Params.filter_min;
    psla = Params.psla;

    pgi_sla_1 = Integrals.IABCsame(:,icomb);  %p(gi | alpha_i = 0, beta_i = 0) (Equation 9)
    pgi_sla_2Amax = 2*Integrals.IA_BC_Amax(:,icomb); pgi_sla_2Amax(pgi_sla_2Amax<0) = 0; %(this is unlikely: should all be > 0
    pgi_sla_2Bmax = 2*Integrals.IA_BC_Bmax(:,icomb); pgi_sla_2Bmax(pgi_sla_2Bmax<0) = 0; % but I did find one floating point error)
    pgi_sla_2Cmax = 2*Integrals.IA_BC_Cmax(:,icomb); pgi_sla_2Cmax(pgi_sla_2Cmax<0) = 0; %p(gi | alpha_i = 1, beta_i = 0) Equation (8)

    pgi_sym = [pgi_sla_2Amax pgi_sla_2Bmax pgi_sla_2Cmax];

    pgi_mas_Tmin = 3*[Integrals.IABC_Amin(:,icomb) Integrals.IABC_Bmin(:,icomb) Integrals.IABC_Cmin(:,icomb)]; %p(gi | beta_i = 1, cell type T is a leaf)
    pgi_mas_T = (1/2)*[pgi_mas_Tmin(:,2)+pgi_mas_Tmin(:,3) pgi_mas_Tmin(:,1)+pgi_mas_Tmin(:,3) pgi_mas_Tmin(:,2)+pgi_mas_Tmin(:,1)]; %p(gi | beta_i = 1, topology T) Equation (4)

    pgi_sla = psla(1)*pgi_sla_1 + psla(2)*pgi_sla_2Amax + psla(3)*pgi_sla_2Bmax + psla(4)*pgi_sla_2Cmax; % p(g_i | beta_i = 0) See Equation (14)

    if filter_min %some additional coarse filtering if option chosen
        [~,imax] = max(data,[],2); 
        expressed = data>log_asym_min;
        sym_expressed = data>log_sym_min;
        pgi_mas_T = pgi_mas_T.*expressed.*((sum(expressed,2)>1)*ones(1,3)); %transition genes above expression threshold in at least 2 cell types
        pgi_sym = pgi_sym.*sym_expressed; %marker genes above expression threshold "marked" cell type
        for i=1:3
            pgi_sym(:,t) = pgi_sym(:,t).*(imax==t); %marker gene most highly expressed in "marked" cell type
        end
    end

    pgi_mas_0 = (1/3)*sum(pgi_mas_T,2); %p(g_i | b_i = 1, T = null) Equation (6), Equation (7)
    pgi_mas_T = [pgi_mas_T pgi_mas_0]; %p(g_i | b_i = 1, T)       


end

