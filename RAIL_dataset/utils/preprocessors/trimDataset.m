function dataset_trimmed = trimDataset(dataset, start_cut, end_cut)


dataset_trimmed = dataset;
for i=1:length(dataset)
    dataset_trimmed(i).time = dataset(i).time((1+start_cut):(end-end_cut));
    dataset_trimmed(i).time = dataset_trimmed(i).time - dataset_trimmed(i).time(1);
    
    
    dataset_trimmed(i).pos =  dataset(i).pos((1+start_cut):(end-end_cut),:);
    
    if isfield(dataset, 'vel')
       dataset_trimmed(i).vel =  dataset(i).vel((1+start_cut):(end-end_cut),:);
    end

    if isfield(dataset, 'acc')
       dataset_trimmed(i).acc =  dataset(i).acc((1+start_cut):(end-end_cut),:);
    end
end

end