% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script reproduces the examples in Section 7.2
% The plots are saved to the current directory
%
% Parameters:
% 
% data_type - the dataset to use
% 'mnist' - Uses the MNIST dataset
% 'poker' - Uses the poker dataset
% 'wine' -  - Uses the wine dataset
% 'super' - Uses the superconductor dataset
% erorr_type - the error metric to use
% 'reconstruction' - || K - K_app || / || K ||
% 'angle' - principal angle between subspaces
% 'subspace' - || UU' - VV' || / || UU'||
% n - the number of samples to use
% n_eigs - the number of eigenpairs to calculate
% automatic_n_eigs - if 1, sets n_eigs to account for the # of componenets that account for 90% of the energy, with a maximum of 10
% n_experiments = # of experiments to perform
% kernel_type - the kernel function
% -1 - the kernel used for the p-band extension with parameter signma
% 0 - Gaussian kernel with parameter sigma
% 1 - linear kernel
% 2 - quadratic kernel
% sigma = parameter for the kernel
% to_gl - if 1, normalizes the kernel to become the symmetric nomralized
% graph Lapacian
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear;

% main parameters

data_type = 'wine';
error_type = 'reconstrucion';
n = 500;
n_eigs = 1;
automatic_n_eigs = 1;
n_experiments = 20;
n_sigmas = 20;
kernel_type = 0;
to_gl = 1;
to_sparse = 0;
percent = 0.2;
big_number_for_stability = 10000000000000;
m = [];

% define sigma_std, the std in which sigma values will be generate to
% produce a range of Hoyer scores
if strcmp(data_type, 'mnist') == 1
    sigma_std = 1000;
elseif strcmp(data_type, 'poker') == 1
    sigma_std = 10;
elseif strcmp(data_type, 'wine') == 1
    sigma_std = 100;
elseif strcmp(data_type, 'super') == 1
    sigma_std = 100;
end

for experiment = 1:n_experiments
    
    fprintf('exp = %d / %d \n', experiment, n_experiments)
    X = load_data(data_type, n, m);
    for i = 1:n_sigmas
        
        % generate Kernels until their Hoyer score falls in the bin, to
        % cover as much hoyer scores between 0 and 1
        bin_start = 0.2 * mod(i,5);
        bin_end = bin_start + 0.2;
        sigma = rand(1) * 10;
        K = gaussian_kernel(kernel_type, X, sigma, to_sparse, to_gl);
        the_hs = hoyer_score(K);
        idx = 0;
        while (the_hs < bin_start || the_hs > bin_end) && idx < 50
            sigma = rand(1) * sigma_std;
            K = gaussian_kernel(kernel_type, X, sigma, to_sparse, to_gl);
            the_hs = hoyer_score(K);
            idx = idx + 1;
        end

        hs(1, experiment, i) = the_hs;
        hs(2, experiment, i) = the_hs;
        hs(3, experiment, i) = the_hs;
        hs(4, experiment, i) = the_hs;
        
        
        if automatic_n_eigs
            the_eigs = sort(real(eig(K)), 'descend');
            n_eigs = min(5,find(cumsum(the_eigs)/sum(the_eigs) > 0.9, 1));
            %fprintf('automatic n eigs = %d\n', n_eigs);
        end
        
        % generate appoximations: l-band, p-band, sparse, etc.

        [Q_real, D_real] = eigs(K, n_eigs, big_number_for_stability);

        % l band
        l = percent * n;
        K_s = K(1:l, 1:l); K_s = [K_s , zeros(l, n - l) ; zeros( n - l, n)];
        [Q_0, D_0] = eigs(K_s, n_eigs, big_number_for_stability); Q_0 = real(Q_0); D_0 = real(D_0);
        dK = K - K_s;
        mu = 0;

        [Q_l_block, D_l_block] = update_pert(Q_0, D_0, dK, mu);

        errors(1,1,experiment, i) = return_error('projection', Q_real, D_real, Q_l_block, D_l_block, K);
        errors(2,1,experiment, i) = return_error('best_app_2_norm', Q_real, D_real, Q_l_block, D_l_block, K);
        errors(3,1,experiment, i) = return_error('angle', Q_real, D_real, Q_l_block, D_l_block, K);
        
        % block-diagonal
        q = floor((l^2 / 2)^0.5);
        mu = 0;
        
        % block 1
        K_s = K(1:q, 1:q); K_s = [K_s , zeros(q, n - q) ; zeros( n - q, n)];
        [Q_0, D_0] = eigs(K_s, n_eigs, big_number_for_stability); 
        dK = K - K_s;
        [Q_l_block_diag1, D_l_block_diag1] = update_pert(Q_0, D_0, dK, mu);
        
        % block 2
        K_s = K(end - q + 1:end, end - q + 1:end); K_s = [zeros(n, n - q) , [zeros( n - q, q) ; K_s] ];
        [Q_0, D_0] = eigs(K_s, n_eigs, big_number_for_stability); 
        dK = K - K_s;
        [Q_l_block_diag2, D_l_block_diag2] = update_pert(Q_0, D_0, dK, mu);

        Q_l_block_diag = [Q_l_block_diag1 , Q_l_block_diag2];
        D_l_block_diag = diag([diag(D_l_block_diag1) ; diag(D_l_block_diag2)])/2;
 
        errors(1, 2, experiment, i) = return_error('kernel_2_norm', Q_real, D_real, Q_l_block_diag, D_l_block_diag, K);
        errors(2, 2, experiment, i) = return_error('best_app_2_norm', Q_real, D_real, Q_l_block_diag, D_l_block_diag, K);
        errors(3, 2, experiment, i) = return_error('angle', Q_real, D_real, Q_l_block_diag, D_l_block_diag, K);


        % p band
        p = floor((((2 *n - 1) - ((2*n - 1)^2 - 4 * (l^2 - n))^0.5)/2)) + 1;
        K_s = make_tridigonal(K, p);
        [Q_0, D_0] = eigs(K_s, n_eigs, big_number_for_stability); 
        dK = K - K_s;
        mu = 0;
        [Q_p_pand, D_p_pand] = update_pert(Q_0, D_0, dK, mu);
        
        errors(1, 3, experiment, i) = return_error('kernel_2_norm', Q_real, D_real, Q_p_pand, D_p_pand, K);
        errors(2, 3, experiment, i) = return_error('best_app_2_norm', Q_real, D_real, Q_p_pand, D_p_pand, K);
        errors(3, 3, experiment, i) = return_error('angle', Q_real, D_real, Q_p_pand, D_p_pand, K);
        
        % sparse
        s = l^2;
        K_s = filter_largest_entries(K, s);
        [Q_0, D_0] = eigs(K_s, n_eigs, big_number_for_stability); 
        dK = K - K_s;
        mu = 0;
        [Q_sparse, D_sparse] = update_pert(Q_0, D_0, dK, mu);

        errors(1,4,experiment, i) = return_error('projection', Q_real, D_real, Q_sparse, D_sparse, K);
        errors(2,4,experiment, i) = return_error('best_app_2_norm', Q_real, D_real, Q_sparse, D_sparse, K);
        errors(3,4,experiment, i) = return_error('angle', Q_real, D_real, Q_sparse, D_sparse, K);
    end
end

% plots
file_name = strcat(data_type,'_sigma_',num2str(sigma));
if strcmp(error_type, 'subspace') == 1
    plot_combined(squeeze(errors(1,:,:,:)), hs, 'subspace', file_name);
elseif strcmp(error_type, 'angle') == 1
    plot_combined(squeeze(errors(3,:,:,:)), hs, 'principal angle', file_name);
else
    plot_combined(squeeze(errors(2,:,:,:)), hs, 'reconstruction', file_name);
end
