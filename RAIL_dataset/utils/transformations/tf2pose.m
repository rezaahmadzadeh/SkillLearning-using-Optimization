function pose = tf2pose(tfMat)
%% mask
R = tfMat(1:3, 1:3);
t = tfMat(1:3, 4);

quat = rotm2quat(R);
trans = t';

%changing to ROS convention
pose = [trans, quat(2:end), quat(1)];

end