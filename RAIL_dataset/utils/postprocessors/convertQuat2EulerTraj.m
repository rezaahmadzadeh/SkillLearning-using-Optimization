function eef_traj_converted = convertEuler2QuatTraj(eef_traj)

eef_traj_converted = zeros(size(eef_traj, 1), 7);
eef_traj_converted(:,1:3) = eef_traj(:,1:3); 

eul = eef_traj(:, 4:6);
quat = eul2quat(eul);

eef_traj_converted(:, 4:7) = [quat(:, 2:end) quat(:, 1)]; % converting to ROS convention 

end