function env = getEnvironments()

%[x,y,z,x,y,z,w]
% NOTE: THESE ARE PEN TIP STATES!
% 
% NOTE: offset of magnet base from pen tip: [0,0.04,0];
penTip2MagOffset = [0,0.04,0,0,0,0,0];

env = containers.Map;

%% REACHING
env('REACHING/init/A') = [0.391, 0.292, 0.239, 0.766, -0.642, -0.000, 0.029];%%[0.518, -0.288, 0.465, 0.751, -0.659, -0.034, 0.006];
env('REACHING/init/B') = [0.468, -0.123, 0.298, 0.727, -0.685, -0.019, -0.026];
env('REACHING/init/C') = [0.742, -0.286, 0.370, 0.877, -0.480, -0.002, -0.009];
env('REACHING/init/D') = [0.595, 0.020, 0.152, 0.760, -0.650, -0.001, 0.018];

% Pen_tip coordinates of EEF at the goal location
env('REACHING/goal/A') = [0.682, -0.028, 0.024, 0.775, -0.632, -0.018, 0.005];
env('REACHING/goal/B') = [0.610, 0.455, 0.021, 0.751, -0.660, -0.002, 0.006];
env('REACHING/goal/C') = [0.585, 0.179, 0.020, 0.743, -0.670, 0.000, 0.010];

env('REACHING/obj/A') = [0.697, -0.032, 0.018, -0.030, 0.011, -0.673, 0.739];
env('REACHING/obj/B') = [0.619, 0.453, 0.020, -0.025, -0.012, -0.675, 0.738];
env('REACHING/obj/C') = [0.592, 0.180, 0.025, -0.039, 0.000, -0.673, 0.738];

%% PUSHING
% env('PUSHING/init/A') = [0.588, -0.302, 0.264, 0.894, -0.446, -0.037, 0.023];
% env('PUSHING/init/B') = [0.424, -0.066, 0.424, 0.784, -0.618, -0.042, 0.035];
% env('PUSHING/init/C') = [0.395, 0.225, 0.192, -0.677, 0.732, 0.056, -0.046];
% env('PUSHING/init/D') = [0.539, 0.091, 0.112, 0.732, -0.678, -0.029, 0.063];
% env('PUSHING/init/E') = [0.433, -0.197, 0.326, 0.777, -0.629, -0.022, 0.032];
% env('PUSHING/init/F') = [0.447, 0.239, 0.177,0.777, -0.629, -0.022, 0.032];
% env('PUSHING/init/G') = [0.259, 0.500, 0.429,0.777, -0.629, -0.022, 0.032];

env('PUSHING/init/A') = [0.343, -0.041, 0.402, 0.784, -0.618, -0.042, 0.035];
env('PUSHING/init/B') = [0.541, -0.278, 0.333,  0.777, -0.629, -0.022, 0.032];%[0.526, -0.371, 0.337, 0.777, -0.629, -0.022, 0.032];
env('PUSHING/init/C') = [0.653, 0.196, 0.049,0.777, -0.629, -0.022, 0.032];
env('PUSHING/init/D') = [0.259, 0.500, 0.429,0.777, -0.629, -0.022, 0.032];

env('PUSHING/goal/A')  = [0.766, 0.110, 0.035, -0.520, 0.853, 0.031, -0.010];
env('PUSHING/goal/B') = [0.761, 0.310, 0.043, 0.984, -0.179, 0.007, -0.007];
env('PUSHING/goal/C') = [0.744, 0.197, 0.044, 0.746, -0.666, -0.012, 0.004];

env('PUSHING/obj/A') = [0.846, 0.064, 0.078, -0.050, -0.045, 0.870, -0.488];
env('PUSHING/obj/B') = [0.814, 0.393, 0.075, 0.002, 0.049, -0.204, 0.978];
env('PUSHING/obj/C') = [0.842, 0.200, 0.065, -0.045, -0.014, -0.682, 0.730];

%% PICKING
env('PICKING/init/A') = [0.518, -0.287, 0.466, 0.751, -0.660, -0.033, 0.005];
env('PICKING/init/B') = [0.468, -0.123, 0.298, 0.727, -0.686, -0.019, -0.026];
env('PICKING/init/C') = [0.742, -0.287, 0.370, 0.877, -0.480, -0.002, -0.009];

% Magnet_tip coordinates of EEF at the goal location
env('PICKING/goal/A') = offsetPoseRef2BasePositionOnly([0.544, 0.214, 0.158, 0.995, 0.081, -0.040, 0.036], penTip2MagOffset);
env('PICKING/goal/B') = offsetPoseRef2BasePositionOnly([0.704, 0.275, 0.154, 0.990, 0.125, -0.006, 0.057], penTip2MagOffset);
env('PICKING/goal/C') = offsetPoseRef2BasePositionOnly([0.600, 0.324, 0.160, 0.998, 0.065, -0.013, -0.018], penTip2MagOffset);

env('PICKING/obj/A') = [0.557, 0.193, 0.128, 0.014, -0.001, -0.666, 0.746];
env('PICKING/obj/B') = [0.714, 0.256, 0.128, -0.010, -0.007, -0.653, 0.757];
env('PICKING/obj/C') = [0.617, 0.305, 0.125, 0.010, -0.040, -0.675, 0.737];
%% PRESSING
% env('PRESSING/init/A') = [0.567, -0.360, 0.470, 0.719, -0.695, -0.004, 0.011];
% env('PRESSING/init/C') = [0.741, 0.030, 0.353, 0.743, -0.669, -0.019, 0.026];

env('PRESSING/init/A') = [0.479, 0.164, 0.476, 0.974, -0.225, -0.038, 0.014];
env('PRESSING/init/B') = [0.637, -0.1177,0.373, 0.743, -0.669, -0.019, 0.026];
env('PRESSING/init/C') = [0.624, -0.113, 0.206, 0.950, -0.313, 0.001, 0.008];
env('PRESSING/init/D') =  [0.447, -0.335, 0.424, 0.996, 0.072, -0.053, 0.022];

% Magnet_tip coordinates of EEF at the goal location
env('PRESSING/goal/A') = offsetPoseRef2BasePositionOnly([0.774, -0.279, 0.055, 0.997, 0.058, -0.009, 0.054], penTip2MagOffset);
env('PRESSING/goal/B') = offsetPoseRef2BasePositionOnly([0.700, -0.213, 0.056, 0.996, 0.089, -0.007, 0.030], penTip2MagOffset);
env('PRESSING/goal/C') = offsetPoseRef2BasePositionOnly([0.761, -0.094, 0.051, 0.949, -0.314, -0.027, -0.013], penTip2MagOffset);


env('PRESSING/obj/A') = [0.799, -0.282, 0.056, 0.018, -0.013, -0.663, 0.748];
env('PRESSING/obj/B') = [0.727,-0.270,0.057, 0.052, -0.020, 0.997, -0.053];
env('PRESSING/obj/C') = [0.797, -0.166, 0.057, 0.021, -0.008, -0.821,  0.570];

% [0.798, -0.404, 0.058, 0.021, -0.007, 0.938, -0.346];
%% WRITING
% env('WRITING/init/A') = [0.474, -0.377, -0.012, 0.899, -0.432, -0.074, 0.007];
% env('WRITING/init/B') = [0.474, -0.402, -0.012, 0.898, -0.435, -0.061, 0.022];
% env('WRITING/init/C') = [0.474, -0.444, -0.013, 0.873, -0.479, -0.043, 0.083];
% env('WRITING/init/D') = [0.474, -0.422, -0.012, 0.865, -0.498, -0.036, 0.049];

env('WRITING/init/A') = [0.459, -0.324, -0.010, 0.762, -0.647, 0.015, 0.009];
env('WRITING/init/B') = [0.476, -0.423, -0.010, 0.777, -0.628, -0.001, -0.041];
env('WRITING/init/C') = [0.544, -0.384, -0.010, 0.787, -0.617, -0.003, -0.012];
env('WRITING/init/D') = [0.495, -0.452, -0.010, 0.815, -0.579, -0.010, 0.008];


env('WRITING/goal/A') = [0, 0, 0, 0.899, -0.432, -0.074, 0.007];
env('WRITING/goal/B') = [0, 0, 0, 0.899, -0.432, -0.074, 0.007];
env('WRITING/goal/C') = [0, 0, 0, 0.899, -0.432, -0.074, 0.007];

env('WRITING/obj/A') = [0.430,-0.555,-0.0220,0,0,0,1];%[0, 0, 0, 0.899, -0.432, -0.074, 0.007];
env('WRITING/obj/B') = [0.430,-0.555,-0.0220,0,0,0,1];
env('WRITING/obj/C') = [0.430,-0.555,-0.0220,0,0,0,1];

end


function poseOffseted =  offsetPoseRef2BasePositionOnly(pose, offset)
poseOffseted(1:3) = transformPtObj2BasePositionOnly(offset, pose);
poseOffseted(4:7) = pose(4:end);
end