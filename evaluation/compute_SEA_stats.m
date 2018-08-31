function [meanSEA, stdSEA, SEA] = compute_SEA_stats(demos,repros)
% function to compute SEA statistics

if iscell(demos) && iscell(repros)
    SEA = zeros(length(demos));
    for i = length(demos)
        SEA(i) = swept_error_area(demos{i},repros{i});
    end
    meanSEA = mean(SEA);
    stdSEA = std(SEA);
else
    SEA = swept_error_area(demos,repros);
    meanSEA = SEA;
    stdSEA = [];
end


