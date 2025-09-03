function Pi_hat_fun = Pi_unpre(A_x, A_w, S_r, E_w, J_r, t_k)

    A_x = double(A_x);
    A_w = double(A_w);
    S_r = double(S_r);
    E_w = double(E_w);
    J_r = double(J_r);
    t_k = double(t_k);

    
   %Pi_hat_fun = @(t) (expm(A_x*(t - t_k)) * E_w + integral(@(tau) expm(A_x*(t - tau)) * A_w * expm(S_r*tau) * J_r, t_k, t, 'ArrayValued', true))*inv(J_r)*expm(-S_r*(t - t_k));
   % Pi_hat_fun = @(t) (expm(A_x*(t - t_k)) * E_w * inv(J_r) + integral(@(tau) expm(A_x*(t - tau)) * A_w * expm(S_r*tau) , t_k, t, 'ArrayValued', true))*expm(-S_r*(t - t_k));

    Pi_hat_fun = @(t) (expm(A_x*(t - t_k)) * E_w * inv(J_r) + integral(@(tau) expm(A_x*(t - tau)) * A_w * expm(S_r*tau) , t_k, t, 'ArrayValued', true));

end
