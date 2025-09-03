function Pi_hat = Pi_unpre_matrix(A_x, A_w, S_r, E_w, J_r, t, t_k)

    A_x = double(A_x);
    A_w = double(A_w);
    S_r = double(S_r);
    E_w = double(E_w);

    A_bar = [A_x, A_w; zeros(2, 1), S_r];
    
    A_bar_exp = expm(A_bar*(t-t_k));

    
    % Pi_hat = A_bar_exp * [E_w; J_r] * inv(J_r) * expm(-S_r*(t-t_k));
    appoggio = A_bar_exp * [E_w; J_r];
    appoggio = appoggio(1:size(A_x, 1), 1:size(E_w, 2));
    Pi_hat = appoggio * inv(J_r) * expm(-S_r*(t-t_k));
end
