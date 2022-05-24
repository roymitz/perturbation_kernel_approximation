clear;
n = 1000;
m = 10; % number of eigenpairs to use

% generate matrix we wish to approximate
Q_true = orth(randn(n,m));
D_true = diag(linspace(1,0.5,m));
A = Q_true * D_true * Q_true';

% sample a submatrix of A - could be *any* symmetric submatrix of A
A0 = A(1:200, 1:200);
A0 = [A0, zeros(200, 800) ; zeros(800,1000)];

% calc eigendecomposition of A0
[Q0,D0] = eigs(A0, m);

% perform extensions
dA = A - A0;
mu = 0;
%mu = (trace(A0) - sum(sum(D0)))/(n-m);
[Q1,D1] = update_pert(Q0, D0, dA, mu);
[Q2,D2] = update_pert_second_order(Q0, D0, A, dA, mu);

% calc errors
err0 = norm(A - A0)/norm(A);
err1 = norm(A - Q1 * D1 * Q1')/norm(A);
err2 = norm(A - Q2 * D2 * Q2')/norm(A);
fprintf('Error of A0: %f\n', err0);
fprintf('Error of first order approximation: %f\n', err1);
fprintf('Error of second order approximation: %f\n', err2);
