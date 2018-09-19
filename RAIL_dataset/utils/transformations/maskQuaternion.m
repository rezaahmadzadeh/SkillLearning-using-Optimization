function quat_masked_ros = maskQuaternion(quat_ros, mask)
% WARNING: Gimbal lock can be a problem. BE CAREFUL!

%NOTE: zero at a location in the mask would negate the rotation in that axis

quat_matlab = [quat_ros(end), quat_ros(1:3)];
eul = quat2eul(quat_matlab);           
eul_masked = fliplr(mask).*eul;  % default is ZYX for euler angles in matlab
quat_masked = eul2quat(eul_masked);
quat_masked_ros = [quat_masked(2:end), quat_masked(1)];

end