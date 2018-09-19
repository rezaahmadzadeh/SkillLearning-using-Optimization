function dataset_monotonic = monotonizeDataset(dataset, tol)

if nargin==1
   tol = 1e-5;
end

sampling_period = dataset(1).time(2) - dataset(1).time(1);

dataset_monotonic = dataset;

for ind1 = 1:size(dataset,2)
    traj = dataset(ind1).pos;
    time = dataset(ind1).time;
    numDims = size(traj,2);
    
    traj_monotonic = traj(1,:);
    time_monotonic = time(1);
    for i=2:size(traj,1)
        if sum(abs(traj(i,:)-traj_monotonic(end,:)) > repmat(tol, 1, numDims) , 2) > 0
            traj_monotonic = [traj_monotonic; traj(i,:)];
            time_monotonic = [time_monotonic; time(i)];
        end
    end   
    
    dataset_monotonic(ind1).time = time_monotonic;
    dataset_monotonic(ind1).pos = traj_monotonic;

end

end
