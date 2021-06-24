function [sparse_K] = sparsify(K, n_entries_to_leave, type)

    l = size(K,1);
    n_entries_to_leave = min(l * l, n_entries_to_leave);
    if type == 1
        
%         [I,J] = find(K~=0); IDX = [I,J]; 
%         IDX = IDX(IDX(:,1) <= IDX(:,2), :);
%         p = size(IDX,1);
%         II = randperm(p); IDX = IDX(II,:); IDX = IDX(1:min(size(IDX,1),n_entries_to_leave),:);
%         
        Nmax = n_entries_to_leave^2; % get Nmax biggest entries
        [ ~, Ind ] = sort(abs(K(:)),'descend');
        
        %max_values = Avec(1:Nmax);
        [ ind_row, ind_col ] = ind2sub(size(K), Ind);
        IDX = [ind_row, ind_col]; 
        %IDX = IDX(IDX(:,1) <= IDX(:,2), :);
        
        IDX = IDX(1:min(size(IDX,1), n_entries_to_leave), :);

        sparse_K = zeros(size(K));

        for i = 1:size(IDX,1)
            sparse_K(IDX(i, 1), IDX(i, 2)) = K(IDX(i, 1), IDX(i, 2));
            sparse_K(IDX(i, 2), IDX(i, 1)) = K(IDX(i, 2), IDX(i, 1));
        end
        
    elseif type == 2
        sparse_K = zeros(size(K));
        for i = 1:l
            for j = 1:l
                if (abs(j - i)) <= 200
                    sparse_K(i,j) = K(i,j);
                end
            end
        end
    end
    

end

