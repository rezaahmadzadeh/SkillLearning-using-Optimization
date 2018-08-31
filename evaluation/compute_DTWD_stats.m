function [meanDTWD,stdDTWD,DTWD] = compute_DTWD_stats(demos,repros)
% function to compute SEA statistics

if iscell(demos) && iscell(repros)
    DTWD = zeros(length(demos));
    for i = length(demos)
        DTWD(i) = DTW_dis(demos{i},repros{i});
    end
    meanDTWD = mean(DTWD);
    stdDTWD = std(DTWD);
else
    DTWD = DTW_dis(demos,repros);
    meanDTWD = DTWD;
    stdDTWD = [];
    
end


