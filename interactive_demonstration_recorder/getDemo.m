clc,clear,close all
% =================================
% Code: Reza Ahmadzadeh - 2012
% =================================
answer = inputdlg('Enter number of demonstrations:','nDemo',1,{'3'});
nDemo = str2double(answer{1});

h = figure('Position', [800, 100, 600, 600]);
D = cell(1,nDemo);
demos = D; % a variable to store demonstrations in the LASA dataset format
dt = 0.01; % sampling interval

for ii=1:nDemo
    title(['Draw Demo no.' num2str(ii) ' out of ' num2str(nDemo)]);
    xlabel('x_1');ylabel('x_2');
    
    A = axis;
    Xlim = xlim;
    Ylim = ylim;
    P = get_pencil_curve(h);
    set(h,'Position', [800, 100, 600, 600]);
    D{ii} = P;
    axis(A);
    xlim(Xlim);
    ylim(Ylim);
    demos{1,ii}.pos = D{1,ii}'; % LASA dataset format
    demos{1,ii}.time = [dt:dt:size(D{1,ii},1)*dt];
end

%% resample
nPoint = 100;
xx = linspace(0,1,nPoint);
Ds = cell(1,nDemo);
for ii=1:nDemo
    tt = linspace(0,1,size(D{ii},1));
    xs = spline(tt,D{ii}(:,1),xx);
    ys = spline(tt,D{ii}(:,2),xx);
    Ds{ii} = [xs.' ys.'];
end

%% smoothing
Dss = cell(1,nDemo);
span = 0.3;
for ii=1:nDemo
    xs = smooth(Ds{ii}(:,1),span,'rloess');
    ys = smooth(Ds{ii}(:,2),span,'rloess');
    Dss{ii} = [xs ys];
    demos{1,ii}.pos = Dss{1,ii}'; % LASA dataset format
    demos{1,ii}.time = [dt:dt:size(Dss{1,ii},1)*dt];
end

gcf;hold off;clf
hold on
for ii=1:nDemo
    plot(Dss{ii}(:,1),Dss{ii}(:,2),'linewidth',2);
end
title('Demonstrations');
xlabel('x_1');ylabel('x_2');
xlim(Xlim);ylim(Ylim);
grid on;
