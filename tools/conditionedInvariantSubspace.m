function S = conditionedInvariantSubspace(A_bar,B_bar,C_bar,n,q)


    imB = subspaceImage(B_bar);
    ker_C_bar =  null(C_bar, 'rational');
    M1 = imB;
    for i = 1:n+q
        inter = intersectionOperator(M1,ker_C_bar);
        image_A = A_bar*inter;
        M2 = [imB, image_A];
        if(rank([M1, M2]) == rank(M1))
            break;
        end
        M1 = M2;
    end

    S = subspaceImage(M2);
end