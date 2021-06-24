function [handle] = plot_error_type(errors, nnzs, color, type)

    n_experimets = size(errors, 1);
    y_interpolate = [];
    for i = 1:n_experimets
        if type == 1
            plot(nnzs(i,:), errors(i,:), string(color));
        elseif type == 2
            the_grid = linspace(0.05,1,size(errors,2));
            for j = 1:size(errors,2)
                x_experiment = nnzs(i,:);
                y_experiment = errors(i,:);
                y_interpolate = [y_interpolate; interp1(x_experiment, y_experiment, the_grid, 'linear', 'extrap')];
            end
            
            
        end
    end
    if type == 2
        y_interpolate_mean = mean(y_interpolate);
        y_interpolate_std = std(y_interpolate);
        
        
        curve1 = y_interpolate_mean + y_interpolate_std;
        curve2 = y_interpolate_mean - y_interpolate_std;
        x2 = [the_grid, fliplr(the_grid)];
        inBetween = [curve1, fliplr(curve2)];
        h = fill(x2, inBetween, string(color));
        set(h,'facealpha',.2, 'linestyle','none')
    	handle = plot(the_grid, y_interpolate_mean, string(color), 'LineWidth',2);
    end

end

