function [dataset_transformed] = transformBase2ObjDataset(dataset, objPose)

dataset_transformed = dataset;

for ii = 1:length(dataset)
    eef_traj = dataset(ii).pos;
    
    switch nargin
        case 1
            % no extra argument, i.e. use the recorded object pose as
            % reference
            ref_pose = dataset(ii).obj(1,:);
            
            % --- Reza ---
            % in this case the object will be at the origin
            dataset_transformed(ii).obj = [0 0 0 ref_pose(4:7)];
            % ---
        case 2
            % convert every thing w.r.t the initial pose of the robot's ee
            if ischar(objPose) && strcmp(objPose, 'init')
                ref_pose = dataset(ii).pos(1,:);
                % --- Reza ---
                dataset_transformed(ii).obj = transformObj2EEinit(dataset(ii).obj(1,:), ref_pose);
                % ---
            else
                % convert everything w.r.t the given pose
                ref_pose = objPose;
                % --- Reza ---
                dataset_transformed(ii).obj = transformObj2EEinit(dataset(ii).obj(1,:), ref_pose);
                % ---                
            end
        otherwise
            error('wrong number of input arguments');
    end
    
    [eef_traj_transformed] = transformTrajBase2Obj(eef_traj, ref_pose);
    
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



function [eefTrajTransformed] = transformTrajBase2Obj(eefTraj, objPose)

eefTrajTransformed = eefTraj;
objTransform = pose2tf(objPose);

for ii = 1:size(eefTraj,1)
    eefPose = eefTraj(ii,:);
    eefTransform = pose2tf(eefPose);
    eefTransformedObj = inv(objTransform)*eefTransform;
    eefTrajTransformed(ii,:) = tf2pose(eefTransformedObj);
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



function tfMat = pose2tf(pose)
% changing to matlab convention [w,x,y,z]
quat = [pose(7), pose(4:6)];
trans = pose(1:3);

R = quat2rotm(quat);
t = trans';

tfMat = [R, t; zeros(1,3), 1];

end

function pose = tf2pose(tfMat)

R = tfMat(1:3, 1:3);
t = tfMat(1:3, 4);

quat = rotm2quat(R);
trans = t';

%changing to ROS convention
pose = [trans, quat(2:end), quat(1)];

end