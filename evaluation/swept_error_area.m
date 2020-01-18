function SEA = swept_error_area(traj_1,traj_2)
% SweptAreaError.m: Compute the swept error area between two trajectories.
% Refer to "Khansar, Billar, Learning COntrol Lyapunov Function to Ensure
% Stability of Dynamical System-based Robot Reaching Motions" for the
% formula (Eq. 14). CAUTION: Only applicable to 2D trajectories

% NOTE: Trajectories should be of the dimension nxT where n is the
% dimension of the signal and T is the number of samples

SEA = 0; % initial error to zero

% sum over the samples
for t = 1:(size(traj_1,2)-1)
    SEA = SEA+polyarea([traj_1(1,t),traj_1(1,t+1),traj_2(1,t+1),traj_2(1,t)],[traj_1(2,t),traj_1(2,t+1),traj_2(2,t+1),traj_2(2,t)]);
end
