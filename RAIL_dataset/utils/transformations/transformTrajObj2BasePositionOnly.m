function eefTrajTransformed = transformTrajObj2BasePositionOnly(eefTraj, objPose)
eefTrajTransformed = eefTraj;

for ii = 1:size(eefTraj,1)
    eefPt = eefTraj(ii,1:3);
    eefTrajTransformed(ii,1:3) = transformPtObj2BasePositionOnly(eefPt, objPose);

end

end
