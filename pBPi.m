function prob = pBPi(x1, alp, bet, lam1, lam2)

    x1 = x1 / lam2;
    
 
    function res = fn3(x1, m)
            res = 1/factorial(x1)*(m).^(x1).*exp(-m);
    end
    
 
    n_points = 10;
    alpha_jacobi = bet - 1;
    beta_jacobi = alp - 1;
    [nodes, weights] = jacobi_gauss(n_points, alpha_jacobi, beta_jacobi);

 
    gs = 0;
    for k = 1:length(nodes)
        t = nodes(k);
        w = weights(k);
        m_val = lam1 * (1 + t) / 2;
        cond_prob = fn3(x1, m_val);
        gs = gs + w * cond_prob;
    end
    
 
    prob = 1 / beta(alp, bet) * 2^(-alp-bet+1) * gs;
    prob = max(0, min(1, prob));  
end