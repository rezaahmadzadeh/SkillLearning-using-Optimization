function dataset_converted = convertQuat2EulerDataset(dataset)

dataset_converted = dataset;

for i = 1:length(dataset)
    dataset_converted(i).pos = zeros(size(dataset(i).pos, 1), 6); 
    quat = [dataset(i).pos(:, end), dataset(i).pos(:, 4:6)]; %converting to matlab convention
    eul = quat2eul(quat); % using the default 'ZYZ' convention
    
    dataset_converted(i).pos(:, 1:3) = dataset(i).pos(:,1:3);
    dataset_converted(i).pos(:, 4:6) = eul;
    
end


end