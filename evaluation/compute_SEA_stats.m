function [meanSEA, stdSEA, SEA] = compute_SEA_stats(demos,repros)
% function to compute SEA statistics

if iscell(demos) && iscell(repros)
    SEA = zeros(1,length(demos));
    for i =1:length(demos)
        SEA(i) = swept_error_area(demos{i},repros{i});
    end
    meanSEA = mean(SEA);
    stdSEA = std(SEA);
else
    SEA = swept_error_area(demos,repros);
    meanSEA = SEA;
    stdSEA = [];
end


