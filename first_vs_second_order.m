clear;
type = 2;
n = 1000;
m = 10;

font_size = 16;
file_name = strcat('pert_theory_', num2str(type));

if type == 1
    e_norm = [0.00001,0.00005,0.0001,0.0005,0.001,0.005, 0.01];
else
    e_norm = [0.002,0.005,0.01,0.02,0.05,0.1,0.2];
end
%A = randn(n); [u,s,v] = svd(A); A = u(:,1:m) * s(1:m, 1:m) * u(:,1:m)' ;
delta = 10;


E = randn(n); [u,s,~] = svd(E); E = u * s * u'; E = E/norm(E);
A = randn(n); [u,s,~] = svds(A, m); A = u * s * u'; A =  A / norm(A);

%E = A; E(1:m,1:m) = 0; A = A - E;


for i = 1:size(e_norm,2)
    
    if type == 1
        A = A + 0.5 * eye(n);
        [Q0,D0] = eigs(A, m);
        mu = (trace(A) - sum(sum(D0)))/(n - m); 
        E_ = e_norm(i) * E;
    else
        [u,s,~] = svd(A);
        s(1:m, 1:m) = diag(sort(1 + rand(m,1),'descend'));
        s(m+1:end, m+1:end) =  diag(e_norm(i) * ones(n-m,1)); 
        A = u * s * u' ;
        E_ = E  * 1e-6;
        [Q0,D0] = eigs(A, m);
        mu = 0;
    end
    
    A_ = A + E_;
    
    [Qreal,Dreal] = eigs(A_, m);

    [Q1, D1] = upadte_pert(Q0, D0, E_, 0, []);
    [Q2, D2] = upadte_pert_second_order(Q0, D0, A, E_, 0, []);
    [Q1_mu, D1_mu] = upadte_pert(Q0, D0, E_, mu, []);
    [Q2_mu, D2_mu] = upadte_pert_second_order(Q0, D0, A, E_, mu, []);
    
    %err(i,1) = (subspace(Qreal,Q1));
    %err(i,2) = (subspace(Qreal,Q2));
    err(i,1) = min(norm(Qreal(:,1) + Q1(:, 1)), norm(Qreal(:,1) - Q1(:, 1)));
    err(i,2) = min(norm(Qreal(:,1) + Q2(:, 1)), norm(Qreal(:,1) - Q2(:, 1)));
    err(i,3) = min(norm(Qreal(:,1) + Q1_mu(:, 1)), norm(Qreal(:,1) - Q1_mu(:, 1)));
    err(i,4) = min(norm(Qreal(:,1) + Q2_mu(:, 1)), norm(Qreal(:,1) - Q2_mu(:, 1)));
end

if type == 1
    figure; hold on; 
    plot(log(e_norm), log(err(:, 1)), 'b');
    plot(log(e_norm), log(err(:, 2)), '--r');
    plot(log(e_norm), log(err(:, 3)), '-+b');
    plot(log(e_norm), log(err(:, 4)), '-+r');
    xlabel('log(||cE||)', 'fontsize', 14);
    ylabel('log(error)', 'fontsize', 14);
    legend({'1st order (mu = 0)','2nd order (mu = 0)','1st order (mu-mean)','2nd order (mu-mean)'}, 'fontsize', font_size, 'Location', 'southeast');
    set(groot,'defaultAxesTickLabelInterpreter','latex'); 
    hold off;
else
    figure; 
    figure('DefaultAxesFontSize',font_size);
    hold on; 
    plot(log(e_norm), log(err(:, 1)), 'b');
    plot(log(e_norm), log(err(:, 2)), '--r');
    xlabel('log(c)' ,'fontsize', font_size);
    set(gca,'fontsize',font_size)
    ylabel('log(error)','fontsize', font_size);
    legend({'1st order','2nd order'}, 'fontsize', font_size, 'Location', 'southeast');
    set(groot,'defaultAxesTickLabelInterpreter','latex'); 
    hold off;
end

saveas(gcf,file_name, 'epsc');