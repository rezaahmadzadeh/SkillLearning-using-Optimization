function eefTrajTransformed = transformTrajBase2ObjPositionOnly(eefTraj, objPose)
eefTrajTransformed = eefTraj;

for ii = 1:size(eefTraj,1)
    eefPt = eefTraj(ii,1:3);
    eefTrajTransformed(ii,1:3) = transformPtBase2ObjPositionOnly(eefPt, objPose);

end

end
