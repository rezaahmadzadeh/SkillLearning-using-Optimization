function [meanHD, stdHD, HD] = compute_Hausdorff_distance_stats(demos,repros)
% function to compute SEA statistics

if iscell(demos) && iscell(repros)
    HD = zeros(1,length(demos));
    for i =1:length(demos)
        HD(i) = HausdorffDist(demos{i}.',repros{i}.');
    end
    meanHD = mean(HD);
    stdHD = std(HD);
else
    HD = HausdorffDist(demos.',repros.');
    meanHD = HD;
    stdHD = [];
end
