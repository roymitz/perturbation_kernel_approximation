function [X] = load_data(data_type, n, m)

    if strcmp(data_type, 'mnist') == 1
        P = randperm(60000);
        X = loadMNISTImages('Datasets/train-images.idx3-ubyte');X = X';X = X(P,:);
    elseif strcmp(data_type, 'poker') == 1
        X = csvread('Datasets/poker-hand-training-true.data');   X = X(randperm(size(X,1)), 1:(end - 1));
    elseif strcmp(data_type, 'wine') == 1
        X = load('Datasets\wine_dataset.mat'); X = X.wine; X = X(randperm(size(X,1)),:);
    elseif strcmp(data_type, 'super') == 1
        X = load('Datasets\superconductor_dataset.mat'); X = table2array(X.X); X = X(randperm(size(X,1)),:);
    else
        X = zeros(n,m);

    end
    
    X = X(1:n, :);
    X = zscore(X);
    
end

