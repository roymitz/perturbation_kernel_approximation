function [ K ] = gaussian_kernel(kernel_type, X, sigma , to_sparse, to_gl)

n = size(X, 1);
K = zeros(n);
for i = 1:n
    for j = i:n
        if kernel_type == -1
            x = 2 + rand();
            if i == j
                K(i,j) = 1 + randn() * 0.00001;
            else
                K(i,j) = 10^(-(abs(i-j) + randn() * 0.00001) / sigma);
                %K(i,j) = 1/abs(i - j)^sigma + randn() * 0.0001;
            end
            %K(i,j) = K(i,j) + randn() * 0.000001;
        elseif kernel_type == 0
            K(i,j) = exp(-norm(X(i,:) - X(j,:))^2 / sigma);
        elseif kernel_type == 1
            K(i,j) = (1 + X(i,:) * X(j,:)')^1;
        elseif kernel_type == 2
            K(i,j) = (1 + X(i,:) * X(j,:)')^2;
        elseif kernel_type == 3
            K(i,j) = (1 + X(i,:) * X(j,:)')^3;
        end
        K(j,i) = K(i,j);
    end
end

if to_sparse
    %K(abs(K) < prctile(abs(K(:)), 95)) = 0;
    K(abs(K) < 1e-10) = 0;
end

if to_gl
    D = diag(sum(K).^-0.5);
    K = D * K * D;
end

%K = K - min(eig(K)) * eye(n);

end

