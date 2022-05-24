function [err] = return_error(err_type, Q_real, D_real, Q_method, D_method, K_real)

    if strcmp(err_type, 'angle') == 1
        err = rad2deg(subspace(Q_real, Q_method));
    elseif strcmp(err_type, 'best_app_2_norm') == 1
        besk_app = Q_real * D_real * Q_real';
        K_method = Q_method * D_method * Q_method';
        err = norm(besk_app - K_method, 2)/norm(besk_app, 2);
    elseif strcmp(err_type, 'kernel_2_norm') == 1
        K_method = Q_method * D_method * Q_method';
        err = norm(K_real - K_method, 2)/norm(K_real, 2);
    elseif strcmp(err_type, 'kernel_F_norm') == 1
        K_method = Q_method * D_method * Q_method';
        err = norm(K_real - K_method, 'fro')/norm(K_real, 'fro');
    elseif strcmp(err_type, 'projection') == 1
        Q_method = orth(Q_method);
        Q_real = orth(Q_real);
        err = norm(Q_method * Q_method' - Q_real * Q_real', 2) /  norm(Q_real * Q_real', 2);
    elseif strcmp(err_type, 'pdist_mean') == 1
        err_p = diag(pdist2(Q_method',Q_real','euclidean'));
        err_m = diag(pdist2(Q_method',-Q_real','euclidean'));
        err = mean(min([err_p';err_m']));
    elseif strcmp(err_type, 'pdist_max') == 1
        err_p = diag(pdist2(Q_method',Q_real','euclidean'));
        err_m = diag(pdist2(Q_method',-Q_real','euclidean'));
        err = max(min([err_p';err_m']));
    elseif strcmp(err_type, 'pdist_test_mean') == 1
        test_real = Q_real(n_train + 1 : end, :);
        test_app = Q_method(n_train + 1 : end, :);
        err_plus = diag(pdist2(test_real',test_app','euclidean'));
        err_min = diag(pdist2(test_real',-test_app','euclidean'));
        err = mean(min([err_plus;err_min]));
    else
        err = 1;
    end
end

