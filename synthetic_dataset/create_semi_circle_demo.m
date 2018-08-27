function [traj, time] = create_semi_circle_demo(rad, numPts, T)

y = linspace(-rad, rad, numPts);
dt = T/numel(y);
for i = 1:numel(y)
    x(i) = -sqrt(rad^2 - y(i)^2);
end

traj = [x;y];

time = [dt:dt:T];