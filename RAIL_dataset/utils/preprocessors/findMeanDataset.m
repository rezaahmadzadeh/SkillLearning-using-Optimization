function datasetMean = findMeanDataset(dataset)

nPts = size(dataset(1).pos,1);
nDims = size(dataset(1).pos,2);
nDemos = length(dataset);

datasetTensor = zeros(nPts, nDims, nDemos);

for ii=1:length(dataset)
    datasetTensor(:,:,ii) = dataset(ii).pos;
end

datasetMean.pos = mean(datasetTensor,3);
datasetMean.time = dataset(1).time;

end