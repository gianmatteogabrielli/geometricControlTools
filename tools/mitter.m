function K = mitter(A, B, old, new)
    m = size(B, 2);
    n = size(B, 1);
    K = zeros(n, m);
    P = ctrb(A, B);
    tollerance = 1e-6;
    if(rank(P) ~= n)
        disp("Non-reachable system.");
        K = [];
    else
        if(length(unique(eig(A))) ~= length(eig(A)))
            disp("There are no n distinct eigenvalues");
            K = [];
            
        else
        %     syms k [1 n]
        %     eqn = k*A-old*k == 0;
        %     sol = solve([eqn, k~= 0], k);
        %     v_a = subs(k, sol);
        %     clear k
        %     syms k [m, 1]
        %     eqn = v_a*B*k == new-old;
        %     sol = solve(eqn, k);
        %     f_a = subs(k, sol);
        %     K = f_a*v_a;
        % 
        % end
           [v, d] = eig(A');
           for i = 1: size(A, 2)
               if d(i, i) - old <= tollerance;
                   v_a = v(:, i)';
                   break;
               end
           end
           if (v_a*B ~= 0)
             syms k [m, 1]
             eqn = v_a*B*k == new-old;
             sol = solve(eqn, k);
             f_a = subs(k, sol);
             K = f_a*v_a;
           end
           
           
        end

end