function dataset_sampled = downsampleDataset(dataset, factor)

fields = fieldnames(dataset);

dataset_sampled = dataset;
for i=1:length(dataset)
    for j=1:length(fields)
        eval(['dataset_sampled(i).', fields{j} '= dataset(i).', fields{j}, '(1:factor:end,:);']);
    end
end


end
