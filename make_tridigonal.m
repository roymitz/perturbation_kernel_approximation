function [ tridiag_A ] = make_tridigonal( A, k )

    tridiag_A = diag(diag(A,0),0);
    
    for i = 1:(k-1)
        tridiag_A = tridiag_A + diag(diag(A,i),i) + diag(diag(A,-i),-i);
    end

end

