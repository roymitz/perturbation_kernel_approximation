function [Q_pert, D_pert] = update_pert_second_order(Q0, D0, A, dA, mu)

% This function performs the second order truncated perturbation approximation of Section 3
% Assume A = A0 + dA, where A0 eigendecompostiton is known, and we wish to apprximate the eigendecomposition of A
% Input:
% Q0 - The eigenvectors of the matrix we wish to update (A0)
% D0 - The eigenvalues of the matrix we wish to update (A0)
% dA - The perturbation matrix
% mu - parameter for the formula, estimates the unkown eigenvalues of A0
%     
% Output:
% Q_pert - the perturbation apprxiation for the eigenvectors
% D_pert - the perturbation apprxiation for the eigenvalues


    n = size(Q0,1); m = size(Q0,2);
    
    Q_pert = zeros(n,m);
    D_pert = zeros(m,m);
    

    for i = 1:m
        
        
        D_pert(i,i) = (D0(i,i) + Q0(:,i)' * dA * Q0(:,i)) ;
        v = Q0(:,i);
        
        for j = 1:m
            if i~=j
                v = v + ((Q0(:,j)' * dA * Q0(:,i))/(D0(i,i) - D0(j,j))) * (Q0(:,j));
            end
            
        end
        
        r = (eye(n) - Q0 * Q0') * dA * Q0(:,i);
        %norm((1/(D0(i,i) - mu)) * (dA * Q0(:,i) - Q0 * Q0' *dA * Q0(:,i)))
        %fprintf('%f, the norm = %f, %f, %f\n', m, norm((eye(n) - Q0 * Q0')), norm(dA), norm(Q0(:,i)));
        %v = v + (1/(D0(i,i) - mu)) * ((eye(n) - Q0 * Q0') * dA * Q0(:,i));
        %v = v + (1/(D0(i,i) - mu)) * r;
        v = v + (1/(D0(i,i) - mu) - mu/((D0(i,i) - mu)^2)) * r + (1/ ((D0(i,i) - mu)^2)) * A * r;
        
        Q_pert(:,i) = v ;
    end
    Q_pert = normc(Q_pert);
end

