function [Q_pert, D_pert] = update_pert(Q0, D0, dA, mu)

% This function performs the first order truncated perturbation approximation of Section 3
% Assume A = A0 + dA, where A0 partial eigendecompostiton is known, and we wish to approximate the eigendecomposition of A
% Input:
% Q0 - The eigenvectors of the matrix we wish to update (A0) (size n * m)
% D0 - The eigenvalues of the matrix we wish to update (A0) (size m * m)
% dA - The perturbation matrix
% mu - parameter for the formula, estimates the unknown eigenvalues of A0
% (typically, either mu = 0 or mu = mean of the unknown eigevlaues of A0
% - mu = (trace(A0) - sum(sum(D0))) / ( n - m) )
%     
% Output:
% Q_pert - the perturbation approxiation for the eigenvectors
% D_pert - the perturbation approxiation for the eigenvalues

    n = size(Q0,1); m = size(Q0,2);

    Q_pert = zeros(n,m);
    D_pert = zeros(m,m);
    

    for i = 1:m

        tmp = dA * Q0(:,i);
        
        D_pert(i,i) = (D0(i,i) + Q0(:,i)' * tmp) ; 

        v = Q0(:,i);
        
        
        for j = 1:m
            if i~=j
                v = v + ((Q0(:,j)' * tmp)/(D0(i,i) - D0(j,j))) * (Q0(:,j));
            end
            
        end
        
        
        r = tmp - Q0 * (Q0' * tmp); 
        %r = (eye(n) - Q0 * Q0') * (dA * Q0(:,i));
        %norm((1/(D0(i,i) - mu)) * (dA * Q0(:,i) - Q0 * Q0' *dA * Q0(:,i)))
        %fprintf('%f, the norm = %f, %f, %f\n', m, norm((eye(n) - Q0 * Q0')), norm(dA), norm(Q0(:,i)));
        %v = v + (1/(D0(i,i) - mu)) * ((eye(n) - Q0 * Q0') * dA * Q0(:,i));
        v = v + (1/(D0(i,i) - mu)) * r;
        %v = v + (1/(D0(i,i) - mu) + mu/((D0(i,i) - mu)^2)) * r - (1/ ((D0(i,i) - mu)^2)) * A * r;
        Q_pert(:,i) = v ; 
    end

    Q_pert = normc(real(Q_pert));
end

