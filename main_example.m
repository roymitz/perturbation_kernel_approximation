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
error_type = 'reconstruction';
n = 500;
n_eigs = 10;
automatic_n_eigs = 1;
n_experiments = 5;
kernel_type = 0;
sigma = 5;
to_gl = 1;

% other parameters

sparsify_kernel = 0;
big_number_for_stability = 100000000;
m = 1;

%params for extension

n_samples_for_graph = 10;
Ls = floor(linspace(5 * n_eigs,n,n_samples_for_graph));
Ps = floor(linspace(2,n,n_samples_for_graph));
Ss = floor(linspace(5,100,n_samples_for_graph));

% errors
errors = zeros(3, 5, n_experiments, n_samples_for_graph);
nnzs = zeros(5, n_experiments, n_samples_for_graph);

hold on;

for experiment = 1:n_experiments
    
    fprintf('experiment # = %d\n', experiment);
    
    % load data, generate kernel, calc the eigenvectors
    X = load_data(data_type, n, m);
    K = gaussian_kernel(kernel_type, X, sigma, sparsify_kernel, to_gl);
    
    if automatic_n_eigs
        the_eigs = sort(real(eig(K)), 'descend');
        n_eigs = min(10,find(cumsum(the_eigs)/sum(the_eigs) > 0.8, 1));
        fprintf('automatic n eigs = %d\n', n_eigs);
    end
    
    [Q_real, D_real] = eigs(K, n_eigs, big_number_for_stability);
    
    %l-block extension extend using increasing l
    fprintf('performing l-block extensions...\n');
    
    for i = 1:size(Ls,2)
        l = Ls(i);
        K_s = K(1:l, 1:l); K_s = [K_s , zeros(l, n - l) ; zeros( n - l, n)];
        [Q_0, D_0] = eigs(K_s, n_eigs, big_number_for_stability); Q_0 = real(Q_0); D_0 = real(D_0);
        dK = K - K_s;
        mu = 0;
        
        [Q_l_block, D_l_block] = upadte_pert(Q_0, D_0, dK, mu);
        
        errors(1, 1, experiment, i) = return_error('kernel_2_norm', Q_real, D_real, Q_l_block, D_l_block, K);
        errors(2, 1, experiment, i) = return_error('best_app_2_norm', Q_real, D_real, Q_l_block, D_l_block, K);
        errors(3, 1, experiment, i) = return_error('angle', Q_real, D_real, Q_l_block, D_l_block, K);
        nnzs(1, experiment, i) = nnz(K_s)/ nnz(K);
    end
    
    
    %mu-shift extension extend using increasing l, with
    %mu_mean
    fprintf('performing mu-shift extensions...\n');
    
    for i = 1:size(Ls,2)
        l = Ls(i);
        K_s = K(1:l, 1:l); K_s = [K_s , zeros(l, n - l) ; zeros( n - l, n)];
        [Q_0, D_0] = eigs(K_s, n_eigs, big_number_for_stability); 
        dK = K - K_s;
        mu = (trace(K_s) - sum(sum(D_0))) / (size(K_s,2) - size(D_0,2));
        
        [Q_mu_shift, D_mu_shift] = upadte_pert(Q_0, D_0, dK, mu);
        
        errors(1, 2, experiment, i) = return_error('kernel_2_norm', Q_real, D_real, Q_mu_shift, D_mu_shift, K);
        errors(2, 2, experiment, i) = return_error('best_app_2_norm', Q_real, D_real, Q_mu_shift, D_mu_shift, K);
        errors(3, 2, experiment, i) = return_error('angle', Q_real, D_real, Q_mu_shift, D_mu_shift, K);
        nnzs(2, experiment, i) = nnz(K_s) / nnz(K);
    end

    %block-diagonal extension extend using increasingly larger blocks
    fprintf('performing block-diagonal extensions...\n');
    
    for i = 1:size(Ls,2)
        l = Ls(i);
        K_s = K(end - l + 1:end, end - l + 1:end); K_s = [zeros(n, n - l) , [zeros( n - l, l) ; K_s] ];
        [Q_0, D_0] = eigs(K_s, n_eigs, big_number_for_stability); 
        dK = K - K_s;
        mu = 0;
        
        [Q_l_block_diag, D_l_block_diag] = upadte_pert(Q_0, D_0, dK, mu);
        Q_l_block_diag = [Q_l_block_diag , Q_l_block];
        D_l_block_diag = diag([diag(D_l_block_diag) ; diag(D_l_block)])/2;
 
        errors(1, 3, experiment, i) = return_error('kernel_2_norm', Q_real, D_real, Q_l_block_diag, D_l_block_diag, K);
        errors(2, 3, experiment, i) = return_error('best_app_2_norm', Q_real, D_real, Q_l_block_diag, D_l_block_diag, K);
        errors(3, 3, experiment, i) = return_error('angle', Q_real, D_real, Q_l_block_diag, D_l_block_diag, K);
        nnzs(3, experiment, i) = nnz(K_s) / nnz(K) + nnzs(1, experiment, i);
    end
    
    %p-band extension extend using increasing p
    fprintf('performing p-band extensions...\n');
    
    for i = 1:size(Ps,2)
        p = Ps(i);
        K_s = make_tridigonal(K, p);
        [Q_0, D_0] = eigs(K_s, n_eigs, big_number_for_stability); 
        dK = K - K_s;
        mu = 0;
        [Q_p_pand, D_p_pand] = upadte_pert(Q_0, D_0, dK, mu);
        
        errors(1, 4, experiment, i) = return_error('kernel_2_norm', Q_real, D_real, Q_p_pand, D_p_pand, K);
        errors(2, 4, experiment, i) = return_error('best_app_2_norm', Q_real, D_real, Q_p_pand, D_p_pand, K);
        errors(3, 4, experiment, i) = return_error('angle', Q_real, D_real, Q_p_pand, D_p_pand, K);
        nnzs(4, experiment, i) = nnz(K_s) / nnz(K);
    end
    
    
    %sparse extension - extend using increasing # samples
    fprintf('performing sparse max-values extensions...\n');
    
    for i = 1:size(Ss,2)
        s = Ss(i);
        K_s = filter_largest_entries(K, floor(nnz(K) * (s / 100)));
        nnz_sparse(i) = nnz(K_s)/ nnz(K);
        [Q_0, D_0] = eigs(K_s, n_eigs, big_number_for_stability); 
        dK = K - K_s;
        mu = 0;
        [Q_sparse, D_sparse] = upadte_pert(Q_0, D_0, dK, mu);
        
        errors(1, 5, experiment, i) = return_error('kernel_2_norm', Q_real, D_real, Q_sparse, D_sparse, K);
        errors(2, 5, experiment, i) = return_error('best_app_2_norm', Q_real, D_real, Q_sparse, D_sparse, K);
        errors(3, 5, experiment, i) = return_error('angle', Q_real, D_real, Q_sparse, D_sparse, K);
        nnzs(5, experiment, i) = nnz(K_s) / nnz(K);
    end  
        
%    fprintf('performing sparse random-values extensions...\n');
    
%     for i = 1:size(Ss,2)
%         s = Ss(i);
%         K_s = filter_random_entries(K, floor(nnz(K) * (s / 100)));
%         nnz_sparse_rand(i) = nnz(K_s)/ nnz(K);
%         [Q_0, D_0] = eigs(K_s, n_eigs, big_number_for_stability); 
%         dK = K - K_s;
%         mu = 0;
%         [Q_sparse, D_sparse] = upadte_pert(Q_0, D_0, dK, mu, n_eigs);
% 
%         err_sparse_rand(i) = return_error(error_type, Q_real, D_real, Q_sparse, D_sparse);
%     end  
    
    
end

% plots
file_name = strcat(data_type,'_sigma_',num2str(sigma));
if strcmp(error_type, 'subspace') == 1
    plot_combined(squeeze(errors(1,:,:,:)), nnzs, 'subspace', file_name);
elseif strcmp(error_type, 'angle') == 1
    plot_combined(squeeze(errors(3,:,:,:)), nnzs, 'principal angle', file_name);
else
    plot_combined(squeeze(errors(2,:,:,:)), nnzs, 'reconstruction', file_name);
end