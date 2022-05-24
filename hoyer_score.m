function [score] = hoyer_score(v)
    vv = v(:);
    N = size(vv,1);
    score = (N^0.5 - norm(vv,1)/norm(vv,2)) * (N^0.5 - 1)^-1;
end

