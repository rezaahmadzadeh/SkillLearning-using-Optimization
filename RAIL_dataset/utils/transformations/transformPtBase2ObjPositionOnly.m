function ptTransformed = transformPtBase2ObjPositionOnly(pt, objPose)
quatFixed = [0,0,0,1];
pose = [pt(1:3), quatFixed];

objTFMat = pose2tf(objPose);
poseTFMat = pose2tf(pose);
poseNewTFMat = objTFMat \ poseTFMat;
poseTransformed = tf2pose(poseNewTFMat);
ptTransformed = poseTransformed(1:3);

end
