function dataset_diffed =  smoothDatasetJA(dataset, num_samples, lambda)

dataset_diffed = dataset;

for jj=1:numel(dataset)
    time = dataset(jj).time;
    pos = dataset(jj).pos;
    
%     M=200; % number of samples per 5 via-points
%     N = round(size(pos,1)*M/5)
%     num_samples = 500;
%     time_diffed = linspace(time(1), time(end), num_samples).';
%     [trj, ~] = min_jerk_todorov(pos(:,1:3), num_samples,  [], [], []);
%   
    time = linspace(time(1), time(end), size(time,1)).';
    time_diffed = linspace(time(1), time(end), num_samples).';
    trj = pos(:,1:3);
    
    [x, y, z]=filter_JA_modified(trj, lambda, time, time_diffed);
        
    dataset_diffed(jj).time = time_diffed;
    dataset_diffed(jj).pos = [x(:,1), y(:,1), z(:,1)];
    dataset_diffed(jj).vel = [x(:,2), y(:,2), z(:,2)];
    dataset_diffed(jj).acc = [x(:,3), y(:,3), z(:,3)];
    
end

end