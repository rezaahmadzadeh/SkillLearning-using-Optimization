function dataset_fixed = fixPtDataset(dataset, fix_idx)

if fix_idx~=Inf
    fixed_pt = dataset(1).pos(fix_idx, :);
    
    dataset_fixed = dataset;
    for i=1:size(dataset,1)
        diff = fixed_pt - dataset(i).pos(fix_idx,:);
        for j = 1:size(dataset(i).pos,1)
            dataset_fixed(i).pos(j,:) = dataset(i).pos(j,:) + diff;
        end
    end
    
else
    dataset_fixed = dataset;
end


end