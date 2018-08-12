function data = restructureData(S, numDemos, numPoints)
%% prepare the data for GMM
data = zeros(3, numDemos*numPoints);
for ii = 1:numDemos
    data(:,1+(ii-1)*numPoints:ii*numPoints) = [S.demos{ii}.t; S.demos{ii}.pos];
end
end