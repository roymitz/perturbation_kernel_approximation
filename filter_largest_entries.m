function [K_s] = filter_largest_entries(K, n_entries)
    K_s = zeros(size(K));
    [~, sortIndex] = sort(K(:), 'descend');  % Sort the values in descending order
    maxIndex = sortIndex(1:n_entries);  % Get a linear index into A of the 5 largest values
    K_s(maxIndex) = K(maxIndex);
end

