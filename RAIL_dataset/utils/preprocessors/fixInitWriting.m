function dataset_trimmed = fixInitWriting(dataset, x_init)


dataset_trimmed = dataset;
for i=1:length(dataset)
    idx = find(dataset(i).pos(:,1)>=x_init);
    
    dataset_trimmed(i).time = dataset(i).time(idx);
    dataset_trimmed(i).time = dataset_trimmed(i).time - dataset_trimmed(i).time(1);
    
    dataset_trimmed(i).pos = dataset_trimmed(i).pos(idx,:);
    
    if isfield(dataset, 'vel')
       dataset_trimmed(i).vel =  dataset(i).vel(idx,:);
    end

    if isfield(dataset, 'acc')
       dataset_trimmed(i).acc =  dataset(i).acc(idx,:);
    end
end

end