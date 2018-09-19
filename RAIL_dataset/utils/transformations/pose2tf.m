function tfMat = pose2tf(pose)

% changing to matlab convention [w,x,y,z]
quat = [pose(7), pose(4:6)];
trans = pose(1:3);

R = quat2rotm(quat);
t = trans';

tfMat = [R, t; zeros(1,3), 1];

end

