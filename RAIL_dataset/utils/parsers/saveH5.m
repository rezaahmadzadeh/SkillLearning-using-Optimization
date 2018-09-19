function saveH5(dataset, filename, overwriteFlag)
% Save the dataset fto an h5 file in the format 'demo{i}/timeDSName, demo{i}/posDSName, ...'
% Input:
%   dataset: a vector of Demo structs (fields: time, pos, vel, acc)
%   filename: h5 filename with full path and extension
%
% MAKE SURE THE DATASET NAMES ARE CONSISTENT WHEN LOADING TOO
% TODO: MAKE THE NAMING DYNAMIC

if nargin ==1
    filename = 'default_dataset.h5';
    disp('Saving as default_dataset.h5')
    overwriteFlag = 1;
elseif nargin == 2
    overwriteFlag = 0;
end

if exist(filename, 'file') == 2
    if overwriteFlag
        disp('Overwriting existing file!')
        delete(filename)
    else
        error('File Exists -  Not allowed to Overwrite')
    end
end


for j = 1:length(dataset)
    ds_names = fieldnames(dataset(j));
    num_ds = length(ds_names);
    for i =1:num_ds
        sz = eval(['size(dataset(j).', ds_names{i},')']);
        contents = eval(['dataset(j).', ds_names{i}]);
        
        if ~isempty(contents)
            h5create(filename,['/demo', num2str(j,'%02.0f'), '/', ds_names{i}], fliplr(sz));
            h5write(filename,['/demo', num2str(j,'%02.0f'), '/', ds_names{i}], contents');
        end
    end
    
end

end
