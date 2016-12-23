function [prob_asym_topol, prob_sym_topol, prob_marker_topol] = find_gene_probabilities(pgi_mas_T, pgi_sla, pgi_sym, tlik, psla, odds0, data, sym_thresh)
    % asym
    odds_asym_topol = odds0*pgi_mas_T(:,tlik)./pgi_sla; % p(a_i = 1 | g_i, T)/p(a_i = 0 | g_i, T)
    prob_asym_topol = odds_asym_topol./(1+odds_asym_topol); % p(a_i = 1 | g_i, T)
    
    % sym
    p_mas = odds0/(1+odds0); % p(a_i=1)
    p_sla = 1/(1+odds0); % p(a_i=0)
    pg_T = p_mas*pgi_mas_T(:,tlik) + p_sla*pgi_sla; % p(g_i | T)
    psym = pgi_sym(:,tlik); % p(g_i | T, gene i symmetric) i.e. p(g_i | T, b_i = 1)
    pnotsym = pg_T-p_sla*psla(tlik+1)*psym; % p(g_i | T, gene i NOT symmetric). psla(tlik+1) is p(mu_A > mu_BC | a_i = 0) (or equiv)
    psym = psym.*(data(:,tlik)>sym_thresh); % set p(g_i | T, gene i symmetric) to be 0 if mean expression < sym_thresh
    prior_sym_odds = (p_sla*psla(tlik+1))/(1-p_sla*psla(tlik+1)); % p(b_i = 1)/p(b_i = 0);
    odds_sym_topol = prior_sym_odds*psym./pnotsym; % p(b_i = 1 | g_i, T)/p(b_i = 0 | g_i, T)
    prob_sym_topol = odds_sym_topol./(odds_sym_topol+1); % p(b_i = 1 | g_i, T)
    
    %markers
    %pg_T as above
    prob_marker_topol = zeros(size(pgi_sym));
    for t=1:3
        psym = pgi_sym(:,t); % p(g_i | T, gene i symmetric) i.e. p(g_i | T, b_i = 1)
        pnotsym = pg_T-p_sla*psla(t+1)*psym; % p(g_i | T, gene i NOT symmetric). psla(tlik+1) is p(mu_A > mu_BC | a_i = 0) (or equiv)
        psym = psym.*(data(:,t)>sym_thresh); % set p(g_i | T, gene i symmetric) to be 0 if mean expression < sym_thresh
        prior_sym_odds = (p_sla*psla(t+1))/(1-p_sla*psla(t+1)); % p(b_i = 1)/p(b_i = 0);
        odds_marker_topol = prior_sym_odds*psym./pnotsym; % p(b_i = 1 | g_i, T)/p(b_i = 0 | g_i, T)
        prob_marker_topol(:,t) = odds_marker_topol./(odds_marker_topol+1); % p(b_i = 1 | g_i, T)
    end
end