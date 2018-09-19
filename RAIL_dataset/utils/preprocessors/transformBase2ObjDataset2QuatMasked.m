function [dataset_transformed] = transformBase2ObjDataset2QuatMasked(dataset, objPose)

dataset_transformed = dataset;
pos_idx = [1,2,3];
quat_idx = [4,5,6,7];
fixed_quat = [0,0,0,1];


for ii = 1:length(dataset)
    eef_traj = dataset(ii).pos;
    
    switch nargin
        case 1
            % no extra argument, i.e. use the recorded object pose as
            % reference
            ref_pose = dataset(ii).obj(1,:);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % VERY IMP: Removing rotations about x and y axis! The table is
            % flat! Some object location erroneously add x-y rotations
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            mask = [0,0,1];
            quat_obj = ref_pose(4:end);
            
            quat_obj_masked = maskQuaternion(quat_obj, mask);
            ref_pose(4:end) = quat_obj_masked;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            dataset_transformed(ii).obj = ref_pose;
            % ---
        case 2
            % convert every thing w.r.t the initial pose of the robot's ee
            if ischar(objPose)
                if strcmp(objPose, 'init')
                    ref_pose = [dataset(ii).pos(1,pos_idx), fixed_quat];
                    dataset_transformed(ii).obj = ref_pose;
                elseif strcmp(objPose, 'objQuatOnly')
                    ref_pose = [0,0,0, dataset(ii).obj(1,quat_idx)];
                else
                    error('Option doesnt exist');
                end
                % ---
            else
                % convert everything w.r.t the given pose
                ref_pose = objPose;
                dataset_transformed(ii).obj = ref_pose;
                % ---                
            end
        otherwise
            error('wrong number of input arguments');
    end
    
    [eef_traj_transformed] = transformTrajBase2ObjPositionOnly(eef_traj, ref_pose);
    
    dataset_transformed(ii).pos = eef_traj_transformed;
    
end

end
