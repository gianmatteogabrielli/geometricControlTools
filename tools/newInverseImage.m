function inverse_I = newInverseImage(A, V)
    inverse_I = [];
    n = size(A, 2); 
    % in questo caso dobbiamo supporre che la matrice A abbia rango pieno.
    
    if(size(A, 2) == size(A,1))
        if(det(A) ~= 0)
            inverse_I = inv(A)*V; 
        else
            for i = 1:n
                appoggio = ['x', num2str(i)];
                var_names(i) = sym(appoggio);
            end
            v = var_names.';
            %keyboard;
            for i = 1:size(V, 2)
                
                eqns = A*v == V(:, i);
                res = solve(eqns, var_names);
                x = fieldnames(res);
                vector = zeros(n, 1);
                for j = 1:size(V, 1)
                    if(isempty(res.(x{1})))
                        vector = [];
                        break;
                    else
                        vector(j) = res.(x{j});%vector(i) = res.(x{i});
                    end
    
                end
    
                if(rank([inverse_I, vector]) ~= rank(inverse_I))
                    inverse_I = [inverse_I, vector];
                end
    
    
            end
      end

        
    else
        % creo un sistema
        
        for i = 1:n
            appoggio = ['x', num2str(i)];
            var_names(i) = sym(appoggio);
        end
        v = var_names.';
        %keyboard;
        for i = 1:size(V, 2)
            
            eqns = A*v == V(:, i);
            res = solve(eqns, var_names);
            x = fieldnames(res);
            vector = zeros(n, 1);
            for j = 1:n
                if(n == 1)
                    if(res == 0)
                        vector = [];
                        break;
                    end
                end
                if(isempty(res.(x{1})))
                    vector = [];
                    break;
                else
                    vector(j) = res.(x{j});
                end

            end

            if(rank([inverse_I, vector]) ~= rank(inverse_I))
                inverse_I = [inverse_I, vector];
            end


        end
        
    end
    

end