function DTWD = DTW_dis(traj1,traj2,w)
% Computes Dynamic Time Wrapping Distance between two trajectories

% NOTE: Trajectories are preferred in of the dimension Txn where n is the
% dimension of the signal and T is the number of samples
if nargin<3
   w = Inf; 
end

if size(traj1,1)<size(traj1,2) % if the user provides nxT trajectory
    traj1 = traj1.';
    traj2 = traj2.';    
end

n = size(traj1,1);
m = size(traj2,1);
w = max([w,abs(n-m)]);

DTW = zeros(n+1,m+1) +Inf;
DTW(1,1) = 0;

for ind_i = 2:1:n+1
   for ind_j = max(2,ind_i - w):1:min(m+1,ind_i+w) 
       cost = norm(traj1(ind_i-1,:)-traj2(ind_j-1,:));
       DTW(ind_i,ind_j) = cost + min([DTW(ind_i-1,ind_j),DTW(ind_i,ind_j-1),DTW(ind_i-1,ind_j-1)]);      
   end
end
DTWD = DTW(n+1,m+1);

end

