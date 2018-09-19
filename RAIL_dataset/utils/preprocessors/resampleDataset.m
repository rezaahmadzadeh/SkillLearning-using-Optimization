function dataset_resampled = resampleDataset(dataset, num_samples, sampling_method)
% sampling method can either be 'pchip' or 'spline'
%

if nargin == 2
    sampling_method = 'pchip';
end

dataset_resampled = dataset;

for i = 1:length(dataset)
    time = dataset(i).time;
    pos = dataset(i).pos;
    
    time_resampled = linspace(time(1), time(end), num_samples).';
    
    for j = 1:size(pos,2)
        pos_resampled(:,j) = eval([sampling_method, '(time, pos(:,j), time_resampled)']);
    end
    
    dataset_resampled(i).time = time_resampled;
    dataset_resampled(i).pos = pos_resampled;
end


end