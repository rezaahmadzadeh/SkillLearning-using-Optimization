function datasetOffseted = offsetDataset(dataset, posOffset, frameName)
if nargin == 2
    frameName = 'pen_tip';
end

posIdx = [1,2,3];
datasetOffseted = dataset;

for ii = 1:length(dataset)
    datasetOffseted(ii).pos = offsetTrajPositionOnly(posOffset, dataset(ii).pos);
%     datasetOffseted(ii).pos(:,posIdx) = dataset(ii).pos(:,posIdx) ...
%         + posOffsetTranformed;
    
    datasetOffseted(ii).frameName = frameName;
end


end

%%%%%%%%%%
%%
%%%%%%%%%%

function offsetedTraj = offsetTrajPositionOnly(offset, trajPose)

for ii = 1:size(trajPose,1)
    eefPose = trajPose(ii,:);
    offsetedTraj(ii,1:3) = transformPtObj2BasePositionOnly(offset, eefPose);
    offsetedTraj(ii,4:7) = trajPose(ii,4:7);
end

end
