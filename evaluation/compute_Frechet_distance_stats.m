function [meanFD, stdFD, FD] = compute_Frechet_distance_stats(demos,repros)
% function to compute SEA statistics

if iscell(demos) && iscell(repros)
    FD = zeros(1,length(demos));
    for i =1:length(demos)
        FD(i) = DiscreteFrechetDist(demos{i}.',repros{i}.');
    end
    meanFD = mean(FD);
    stdFD = std(FD);
else
    FD = HausdorffDist(demos.',repros.');
    meanFD = FD;
    stdFD = [];
end
