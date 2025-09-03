function [Pi_x, Pi_w, A_x, A_w, E_x, E_w, Omega, S_r, J_r, Gamma_x, Gamma_w, Lambda_x, Lambda_w] = hybridFrancisSolutions(A, B, E, F, P, R, S, J, K, K_d, V, control)
    
    dimV = rank(V);
    n = size(A, 2);
    q = size(S, 2);
    mc = size(K, 1);
    md = size(K_d, 1);
    if(dimV == q)
        disp("--------------------------------------------------------------------------");
        disp("The controlled invariant subspace matches the dimension of the exosystem.");
        disp("--------------------------------------------------------------------------");
        Pi_x = [];
        Pi_w = V(1:n, :);
        S_r = S;
        J_r = J;
        Omega = eye(q);
        
        Gamma_x = K(:, 1:n);
        Gamma_w = K(:, n+1:end);
        Lambda_x = K_d(:, 1:n);
        Lambda_w = K_d(:, n+1:end);
        
        syms sw [n q]
        
        eqns_flow = Pi_w*S_r == A*Pi_w + B*Gamma_x*Pi_w + B*Gamma_w*Omega + P*Omega;
        sol_flow = solve(eqns_flow, sw);

        eqns_jump = Pi_x*sw + Pi_w*J_r == E*Pi_w + F*Lambda_x*Pi_w + F*Lambda_w*Omega + R*Omega;
        sol_jump = solve(eqns_jump, sw);

    else
        
        Pi_x = V(1:n,1:dimV - q);
        Pi_w = V(1:n, dimV - q +1:end);
        S_r = S;
        J_r = J;
        Omega = eye(2);
        
    end
    if(control == 1)
        Gamma_x = K(:, 1:n);
        Gamma_w = K(:, n+1:end);
        Lambda_x = K_d(:, 1:n);
        Lambda_w = K_d(:, n+1:end);
        
        syms sx [size(Pi_x, 2) size(Pi_x, 2)]
        syms sw [n+q-size(Omega, 2) size(Pi_w, 2)]
        
        eqns1 = Pi_x*sx == A*Pi_x + B*Gamma_x*Pi_x;
        eqns3 = Pi_x*sw + Pi_w*S_r == A*Pi_w + B*Gamma_x*Pi_w + B*Gamma_w*Omega + P*Omega;
        vars = [sx(:); sw(:)];
        sol1 = solve([eqns1, eqns3], vars);

        eqns2 = Pi_x*sx == E*Pi_x + F*Lambda_x*Pi_x;
        eqns4 = Pi_x*sw + Pi_w*J_r == E*Pi_w + F*Lambda_x*Pi_w + F*Lambda_w*Omega + R*Omega;
        vars = [sx(:); sw(:)];
        sol2 = solve([eqns2, eqns4], vars);
        [ni, mu] = size(sw);

        A_x = zeros(ni, ni);  
        A_w = zeros(ni, mu);
        E_x = A_x;
        E_w = A_w;
        
        if(ni == 1)
            A_x = sol1.sx1;
            E_x = sol2.sx1;
        else
        % A_x, E_x matrices
            for i = 1:ni
                for j = 1:ni
                    field_name_x = sprintf('sx%d_%d', i, j); % Crea il nome del campo sx{i,j}
                    if isfield(sol1, field_name_x) % Controlla se esiste nel struct
                        A_x(i, j) = sol1.(field_name_x); % Assegna il valore
                    end
                    if isfield(sol2, field_name_x) % Controlla se esiste nel struct
                        E_x(i, j) = sol2.(field_name_x); % Assegna il valore
                    end
                end
            end
        end
        if(mu == 1)
            A_w = sol1.sw1;
            E_w = sol2.sw2;
        else
            % A_w, E_w matrices
            for i = 1:ni
                for j = 1:mu
                    field_name_w = sprintf('sw%d', j); % Crea il nome del campo sx{i,j}
                    if isfield(sol1, field_name_w) % Controlla se esiste nel struct
                        A_w(i, j) = sol1.(field_name_w); % Assegna il valore
                    end
                    if isfield(sol2, field_name_w) % Controlla se esiste nel struct
                        E_w(i, j) = sol2.(field_name_w); % Assegna il valore
                    end
                end
            end
        end

        
        
    else
        syms Gamma_x_v [size(B, 2) size(Pi_x, 1)] % le colonne sono da studiare meglio
        syms Gamma_w_v [size(B, 2) size(Omega, 1)]
        syms Lambda_x_v [size(F, 2) size(Pi_x, 1)]
        syms Lambda_w_v [size(F, 2) size(Omega, 1)]
        syms sx [size(Pi_x, 2) size(Pi_x, 2)] 
        syms sw [size(Pi_x, 2) size(Pi_w, 2)]
        
        eqns1 = Pi_x*sx == A*Pi_x + B*Gamma_x_v*Pi_x;
        eqns3 = Pi_x*sw + Pi_w*S_r == A*Pi_w + B*Gamma_x_v*Pi_w + B*Gamma_w_v*Omega + P*Omega;
        vars = [sx(:); sw(:); Gamma_x_v(:); Gamma_w_v(:)];
        sol1 = solve([eqns1, eqns3], vars);

        % A_x = double([sol1.sx1_1, sol1.sx1_2;sol1.sx2_1,sol1.sx2_2]);
        % A_w = double([sol1.sw1_1, sol1.sw1_2;sol1.sw2_1,sol1.sw2_2]);
        % S_r = S;
        % Gamma_x = double([sol1.Gamma_x_v1_1, sol1.Gamma_x_v1_2, sol1.Gamma_x_v1_3;
        %             sol1.Gamma_x_v2_1, sol1.Gamma_x_v2_2, sol1.Gamma_x_v2_3]);
        % Gamma_w = double([sol1.Gamma_w_v1_1, sol1.Gamma_w_v1_2;
        %             sol1.Gamma_w_v2_1, sol1.Gamma_w_v2_2]);
        % 
  
        eqns2 = Pi_x*sx == E*Pi_x + F*Lambda_x_v*Pi_x;
        eqns4 = Pi_x*sw + Pi_w*J_r == E*Pi_w + F*Lambda_x_v*Pi_w + F*Lambda_w_v*Omega + R*Omega;
        vars = [sx(:); sw(:); Lambda_x_v(:); Lambda_w_v(:)];
        % sol_Lambda_x_v = solve(eqns2(1, :), Lambda_x_v);
        % sol_Lambda_w_v = solve(eqns4(1, :), Lambda_w_v);
        % sol_sx = solve(eqns2(2:end, :), sx);
        % sol_sw = solve(eqns4(2:end, :), sw);
        
        % E_x = double(subs(sx,sol_sx));
        % E_w = double(subs(sw,sol_sw));
        % J_r = J;
        % Lambda_x = double(subs(Lambda_x_v, sol_Lambda_x_v));
        % Lambda_w = subs(Lambda_w_v, sol_Lambda_w_v);
        % Lambda_w = [0, 1];
        sol2 = solve([eqns2, eqns4], vars);
        [ni, mu] = size(sw);

        A_x = zeros(ni, ni);  
        A_w = zeros(ni, mu);
        E_x = A_x;
        E_w = A_w;
        
        if(ni == 1)
            A_x = sol1.sx1;
            E_x = sol2.sx1;
        else
        % A_x, E_x matrices
            for i = 1:ni
                for j = 1:ni
                    field_name_x = sprintf('sx%d_%d', i, j); % Crea il nome del campo sx{i,j}
                    if isfield(sol1, field_name_x) % Controlla se esiste nel struct
                        A_x(i, j) = sol1.(field_name_x); % Assegna il valore
                    end
                    if isfield(sol2, field_name_x) % Controlla se esiste nel struct
                        E_x(i, j) = sol2.(field_name_x); % Assegna il valore
                    end
                end
            end
        end
        if(mu == 1)
            A_w = sol1.sw1;
            E_w = sol2.sw2;
        else
            % A_w, E_w matrices
            for i = 1:ni
                for j = 1:mu
                    field_name_w = sprintf('sw%d_%d', i, j); % Crea il nome del campo sx{i,j}
                    if isfield(sol1, field_name_w) % Controlla se esiste nel struct
                        A_w(i, j) = sol1.(field_name_w); % Assegna il valore
                    end
                    if isfield(sol2, field_name_w) % Controlla se esiste nel struct
                        E_w(i, j) = sol2.(field_name_w); % Assegna il valore
                    end
                end
            end
        end
        
        

        Gamma_x = zeros(size(Gamma_x_v, 1), size(Gamma_x_v, 2));
        Lambda_x = zeros(size(Lambda_x_v, 1), size(Lambda_x_v, 2));
        Gamma_w = zeros(size(Gamma_w_v, 1), size(Gamma_w_v, 2));
        Lambda_w = zeros(size(Lambda_w_v, 1), size(Lambda_w_v, 2));


        for i = 1:size(Gamma_x_v, 1)
            for j = 1:size(Gamma_x_v, 2)
                field_name_x = sprintf('Gamma_x_v%d_%d', i, j); % Crea il nome del campo sx{i,j}
                if isfield(sol1, field_name_x) % Controlla se esiste nel struct
                    Gamma_x(i, j) = sol1.(field_name_x); % Assegna il valore
                end
            end
        end

        for i = 1:size(Gamma_w_v, 1)
            for j = 1:size(Gamma_w_v, 2)
                field_name_w = sprintf('Gamma_w_v%d_%d', i, j); % Crea il nome del campo sx{i,j}
                if isfield(sol1, field_name_w) % Controlla se esiste nel struct
                    Gamma_w(i, j) = sol1.(field_name_w); % Assegna il valore
                end
            end
        end

        for i = 1:size(Lambda_x_v, 1)
            for j = 1:size(Lambda_x_v, 2)
                field_name_x = sprintf('Lambda_x_v%d_%d', i, j); % Crea il nome del campo sx{i,j}
                if isfield(sol2, field_name_x) % Controlla se esiste nel struct
                    Lambda_x(i, j) = sol2.(field_name_x); % Assegna il valore
                end
            end
        end

        for i = 1:size(Gamma_w_v, 1)
            for j = 1:size(Gamma_w_v, 2)
                field_name_w = sprintf('Lambda_w_v%d_%d', i, j); % Crea il nome del campo sx{i,j}
                if isfield(sol2, field_name_w) % Controlla se esiste nel struct
                    Lambda_w(i, j) = sol2.(field_name_w); % Assegna il valore
                end
            end
        end






        

end