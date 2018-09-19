function [dataset_transformed] = transformBase2ObjDataset2(dataset, objPose)

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
            
            % --- Reza ---
            % in this case the object will be at the origin
%             dataset_transformed(ii).obj = [0, 0, 0, fixed_quat];
            dataset_transformed(ii).obj = ref_pose;
            % ---
        case 2
            % convert every thing w.r.t the initial pose of the robot's ee
            if ischar(objPose)
                if strcmp(objPose, 'init')
                    ref_pose = [dataset(ii).pos(1,pos_idx), fixed_quat];
                    % --- Reza ---
                    %                 dataset_transformed(ii).obj = transformObj2EEinit(dataset(ii).obj(1,:), ref_pose);
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
                % --- Reza ---           
%                 dataset_transformed(ii).obj = transformObj2EEinit(dataset(ii).obj(1,:), ref_pose);
                dataset_transformed(ii).obj = ref_pose;
                % ---                
            end
        otherwise
            error('wrong number of input arguments');
    end
    
    [eef_traj_transformed] = transformTrajBase2Obj2(eef_traj, ref_pose);
    
    dataset_transformed(ii).pos = eef_traj_transformed;
    
end

end


function [obj_pose_ee_init] = transformObj2EEinit(obj_pose, ee_init_pose)
% converting the object pose to the initial ee pose
ee_init_poseT = pose2tf(ee_init_pose);
obj_pose_T = pose2tf(obj_pose);
obj_pose_T_ee_init = inv(ee_init_poseT)*obj_pose_T;
obj_pose_ee_init = tf2pose(obj_pose_T_ee_init);
end



function [eefTrajTransformed] = transformTrajBase2Obj2(eefTraj, objPose)

eefTrajTransformed = eefTraj;
objTransform = pose2tf(objPose);
fixed_quat = [0,0,0,1];
pos_idx = [1,2,3];


for ii = 1:size(eefTraj,1)
    eefPose = [eefTraj(ii,:), fixed_quat];
    eefTransform = pose2tf(eefPose);
    eefTransformedObj = inv(objTransform)*eefTransform;
    temp = tf2pose(eefTransformedObj);
    eefTrajTransformed(ii,:) = temp(:,pos_idx);
end
end


function eefTrajTransformed = transformTrajObj2Base(eefTraj, objPose)
eefTrajTransformed = eefTraj;
objTransform = pose2tf(objPose);

for ii = 1:size(eefTraj,1)
    eefPose = eefTraj(ii,:);
    eefTransform = pose2tf(eefPose);
    eefTransformedObj = objTransform*eefTransform;
    eefTrajTransformed(ii,:) = tf2pose(eefTransformedObj);
end

end