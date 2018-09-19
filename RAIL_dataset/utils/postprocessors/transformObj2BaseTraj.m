function eefTrajTransformed = transformObj2BaseTraj(eefTraj, objPose)
eefTrajTransformed = eefTraj;
objTransform = pose2tf(objPose);

for ii = 1:size(eefTraj,1)
    eefPose = eefTraj(ii,:);
    eefTransform = pose2tf(eefPose);
    eefTransformedObj = objTransform*eefTransform;
    eefTrajTransformed(ii,:) = tf2pose(eefTransformedObj);
end

end



function tfMat = pose2tf(pose)
% changing to matlab convention [w,x,y,z]
quat = [pose(7), pose(4:6)];
trans = pose(1:3);

R = quat2rotm(quat);
t = trans';

tfMat = [R, t; zeros(1,3), 1];

end

function pose = tf2pose(tfMat)

R = tfMat(1:3, 1:3);
t = tfMat(1:3, 4);

quat = rotm2quat(R);
trans = t';

%changing to ROS convention
pose = [trans, quat(2:end), quat(1)];

end