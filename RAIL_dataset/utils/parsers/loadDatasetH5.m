function dataset = loadDatasetH5(dataset_dir)
% Loads the dataset from an h5 file in the format 'demo{i}/timeDSName, demo{i}/posDSName, ...'
%   filename: h5 filename with full path and extension
%   ds_name(optional): cell array of strings containing the field names to
%   be loaded. 
% Output:
%   dataset: a vector of Demo structs (fields: time, pos, vel, acc)
%
% MAKE SURE THE DATASET NAMES ARE CONSISTENT WITH THE FILE LOADED

s = dir([dataset_dir, '/', '*.h5']);
file_list = {s.name};

%% load dataset
dataset = [];

for i=1:size(file_list,2)
    Demo = loadDemoH5([dataset_dir, '/', file_list{i}]);
    dataset = [dataset; Demo];
end

end


function Demo = loadDemoH5(filename)

file_info = h5info(filename);
root_ds_name =  file_info.Groups.Name;

eef_ds_prefix = '_robot_limb_right_right_gripper_tip_state';
%eef_ds_prefix = '_robot_limb_right_endpoint_state';

eef_ds_names = {[eef_ds_prefix,'/pose/position/x'], ...
    [eef_ds_prefix,'/pose/position/y'], ...
    [eef_ds_prefix,'/pose/position/z'], ...
    [eef_ds_prefix,'/pose/orientation/x'],...
    [eef_ds_prefix,'/pose/orientation/y'], ...
    [eef_ds_prefix,'/pose/orientation/z'], ...
    [eef_ds_prefix,'/pose/orientation/w'] };


obj_ds_names = {'_obj_pos_array/poses/poses0/position/x', ...
    '_obj_pos_array/poses/poses0/position/y', ...
    '_obj_pos_array/poses/poses0/position/z', ...
    '_obj_pos_array/poses/poses0/orientation/x', ...
    '_obj_pos_array/poses/poses0/orientation/y', ...
    '_obj_pos_array/poses/poses0/orientation/z', ...
    '_obj_pos_array/poses/poses0/orientation/w'};

time_ds_name = {[eef_ds_prefix, '/header/stamp']};

time = extract_dataset(filename, time_ds_name, root_ds_name);
time = time - time(1);

eef_pose = extract_dataset(filename, eef_ds_names, root_ds_name);
obj_pose = extract_dataset(filename, obj_ds_names, root_ds_name);

% time = 0.01*(0:(length(eef_pose)-1)).';

Demo = struct();

Demo.time = time;
Demo.pos = eef_pose;
Demo.obj = obj_pose;

end


function data = extract_dataset(filename, ds_names_cell, root_ds_name)

data = cell(size(ds_names_cell));
for i =1:length(ds_names_cell)
    data{i} = cell2mat(h5read(filename, [root_ds_name,'/', ds_names_cell{i}]));
end

data = cell2mat(data);
end