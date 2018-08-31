function [meanSSE, stdSSE, SSE] = compute_SSE_stats(demos,repros)
% function to compute SEA statistics

if iscell(demos) && iscell(repros)
    SSE = zeros(length(demos));
    for i = length(demos)
        SSE(i) = (sum(sum((demos{i} - repros{i}).^2)));
    end
    meanSSE = mean(SSE);
    stdSSE = std(SSE);
else
    SSE = (sum(sum((demos - repros).^2)));
    meanSSE = SSE;
    stdSSE = [];
end


