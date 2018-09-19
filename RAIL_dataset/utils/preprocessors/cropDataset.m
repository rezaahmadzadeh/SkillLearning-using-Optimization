function dataset_cropped = cropDataset(dataset, tol)

if nargin ==1
   tol = 1e-5;
end

dataset_cropped = dataset;

for ind1 = 1:length(dataset)
    traj = dataset(ind1).pos;
    time = dataset(ind1).time;
   
    numPts = size(traj,1);
    numDims = size(traj,2);
    
    indMin = find(sum(abs(traj - repmat(traj(1,:), numPts, 1)) > repmat(tol, numPts, numDims), 2) > 0, 1, 'first');  
    
    init = traj(1,:);
    initTime = time(indMin-1);
    
    traj = traj((indMin):end, :);
    time = time((indMin):end);
    
    traj = [init; traj];
    time = [initTime; time];
       
    numPts = size(traj,1);
    
    indMax = find(sum(abs(traj - repmat(traj(end, :), numPts, 1)) > repmat(tol, numPts, numDims), 2) > 0, 1, 'last');
    
    goal = traj(end,:);
    goalTime = time(indMax+1);
    
    traj = traj(1:(indMax), :);
    time = time(1:(indMax));
    
    traj = [traj; goal];
    time = [time; goalTime];
      
    time = time - time(1);
    
    dataset_cropped(ind1).time = time;
    dataset_cropped(ind1).pos = traj;
end

end
