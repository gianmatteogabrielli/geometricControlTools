function c_inv = isControlledInvariant(V, A, B, E, F)
    
    a1 = A*V;
    b1 = [V, B];
    if(rank([a1, b1]) ~= rank(b1))
        c_inv = -1;
    else
        e1 = E*V;
        f1 = [F, V];
        if(rank([e1, f1]) ~= rank(f1))
            c_inv = -1;
        else
            c_inv = 0;
        end
    end
end