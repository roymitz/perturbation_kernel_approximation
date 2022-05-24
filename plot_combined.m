function [] = plot_combined(errors, nnzs, title, file_name)
    type = 2;
    font_size = 16;
    n_error_type = size(errors,1);
    colors = {'blue', 'red', 'black', 'green', 'cyan'};
    %f = figure();
    figure('DefaultAxesFontSize',font_size)
    %set(groot,'defaultAxesTickLabelInterpreter','latex');  
    xlim([0,1]);
    ylim([0,1.5])
    xlabel('Hoyer score', 'fontsize', font_size);
    set(gca,'fontsize',font_size)
    err_string = ['error [', title, ']'];
    
    ylabel(err_string, 'fontsize', font_size);
    hold on;
    h = [];
    for i = 1:n_error_type
        if type == 1
            plot_error_type(squeeze(errors(i,:,:)), squeeze(nnzs(i,:,:)), colors(i), 1);
        elseif type == 2
            h(i) = plot_error_type(squeeze(errors(i,:,:)), squeeze(nnzs(i,:,:)), colors(i), 2);
        elseif type == 3
            plot_error_type(squeeze(errors(i,:,:)), squeeze(nnzs(i,:,:)), colors(i), 1);
        end
    end
    hold off;
    if type == 2
        legend(h, {'l-block','block-diagnoal','p-band','sparse'}, 'fontsize', font_size);
        %legend(h, {'l-block','sparse'}, 'fontsize', font_size);
    end

    set(groot,'defaultAxesTickLabelInterpreter','latex'); 
    saveas(gcf,file_name, 'epsc');
end

