function dataset_smooth = smoothDataset(dataset, smoothing_method, smoothing_param)
% check the input options of smooth function
%
%

if nargin == 1
    smoothing_method = 'moving';
    smoothing_param = 10;
end

dataset_smooth = dataset;

for i = 1:length(dataset)
    time = dataset(i).time;
    pos = dataset(i).pos;
    
    time_smooth = time;
    pos_smooth = [];
    for j = 1:size(pos,2)
        pos_smooth(:,j) = smooth(time, pos(:,j), smoothing_param, smoothing_method);
    end
    
    dataset_smooth(i).time = time_smooth;
    dataset_smooth(i).pos = pos_smooth;
end


end