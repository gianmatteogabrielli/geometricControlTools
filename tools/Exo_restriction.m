function N_W = Exo_restriction(N_star, n, q)
        
    X = [eye(n); zeros(q, n)];  
    N_star_X = intersectionOperator(N_star, X);
    N_star_X_perp = null(N_star_X', 'rational'); 
    N_W = intersectionOperator(N_star, N_star_X_perp);
end
