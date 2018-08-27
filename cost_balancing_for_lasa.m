%% Cost-balancing for LASA dataset
% This file trains GMM and GMM-delta on lasa dataset
% computes the best combination of shape and position objectives
% using optimzization


clc, clear, close all

%--------------------------------
% INITIALIZATION
%--------------------------------
% add paths
addpath('LASA_dataset');
addpath('synthetic_dataset');
addpath('encoder');
addpath('encoder/GMM_GMR');
addpath('encoder/GMM_GMR/lib');
addpath('meta_optimization');
addpath('interactive_demonstration_recorder');
% setup CVX toolbox
% run('C:\Users\Reza\Documents\MATLAB\cvx\cvx_setup.m')

% get all skills from the dataset folder
foldername = 'synthetic_dataset'; % folder name containing demos
ext = 'mat';
[skills,nskills] = getAllMatFiles(foldername, ext); % skills{1}

%--------------------------------
% SKILL LOOP (for more than 1 dataset)
%--------------------------------
doDownSampling = 0; % 1 or 0
doTimeAlignment = 1; % 1 or 0
doSmoothing = 1; % 1 or 0
fixedWeight = 1;        %1e9 weight should not be used because the constraint is included in the optimization;

kk = 2; % skill number (choose between 1:nskills)
load(skills{kk});% loads a dataset including a demo cell and an average dt each demo includes pos, t, vel, acc, dt
nbDemos = size(demos,2);            % number of demos
nbStatesPos = 5;                    % number of Gaussian Components (for position)
nbStatesDelta = 5;                  % number of Gaussian Components (for laplacian)
nbNodes = size(demos{1}.pos,2);     % number of points in each demonstrations
nbDims   = size(demos{1}.pos,1);    % number of dimension (2D / 3D)

%--------------------------------
% Time align the demonstrations
%--------------------------------

if doTimeAlignment
    demos = alignDataset(demos,1);
end


%--------------------------------
% DownSample
%--------------------------------
if doDownSampling
    Demos = cell(1,nbDemos);
    stp = floor(nbNodes / floor(nbNodes * 0.10));
    for ii = 1:nbDemos
        Demos{ii} = demos{ii}.pos(:,1:stp:end);
    end
    nbNodes = size(Demos{1},2);
    
else
    for ii=1:nbDemos
        Demos{ii} = demos{ii}.pos;
    end
end

clear demos dt ext foldername stp ii doDownSampling

%--------------------------------
% Smooth the demonstrations
%--------------------------------
if doSmoothing
    for ii=1:nbDemos
        for j = 1:size(Demos{ii},1)
            Demos{ii}(j,:) = smooth(Demos{ii}(j,:));
        end
    end
end

%--------------------------------
% GMM/GMR - in position space for different nbstates for comp reasons
%--------------------------------
Gmms = cell(1,4);                       % to save GMM/GMR results
D1 = zeros(nbDims+1, nbDemos*nbNodes);  % restructuring the data
t = 1:nbNodes;                          % index
D1(1,:) = repmat(t, 1, nbDemos);
for ii=1:nbDemos
    D1(2:3, (ii-1)*nbNodes+1:ii*nbNodes) = Demos{ii};
end

for ns = 4:7
    M = encodeGMM(D1, nbNodes, ns);
    [repro1, expSigma1] = reproGMM(M);
    M.repro = repro1;
    M.expSigma = expSigma1;
    Gmms{1,ns-3} = M;
end
clear D1 t ns M repro1 expSigma1

%--------------------------------
% GMM/GMR - in position space
%--------------------------------
[Mu_x, R_Sigma_x] = trainGMM(Demos, nbDims, nbDemos, nbNodes, nbStatesPos);
% figure;hold on;
% title('GMM');
% for ii=1:nbDemos
%     plot(Demos{ii}(1,:),Demos{ii}(2,:),'color',[0.5 0.5 0.5]);
% end
% plot(Mu_x(1,:),Mu_x(2,:),'r','linewidth',2)
% bound_x = abs(max(Mu_x(1,:)) - min(Mu_x(1,:)))*0.1;
% bound_y = abs(max(Mu_x(2,:)) - min(Mu_x(2,:)))*0.1;
% axis([min(Mu_x(1,:))-bound_x max(Mu_x(1,:))+bound_x min(Mu_x(2,:))-bound_y max(Mu_x(2,:))+bound_y]);
% xticklabels([]);
% yticklabels([]);
% box on; grid on;
% ylabel('x_2','fontname','Times','fontsize',14);
% xlabel('x_1','fontname','Times','fontsize',14);

% clear bound_x bound_y ii

%--------------------------------
% GMM/GMR - in Gradient space
%--------------------------------
[Mu_g, R_Sigma_g, G] = trainGMMG(Demos, nbDims, nbDemos, nbNodes, nbStatesDelta);

%--------------------------------
% GMM/GMR - in Laplace space
%--------------------------------
[Mu_d, R_Sigma_d, L] = trainGMML(Demos, nbDims, nbDemos, nbNodes, nbStatesDelta);

%--------------------------------
% META-OPTIMIZATION
%--------------------------------
% initialization
M.nbDims = nbDims;
M.nbNodes = nbNodes;
M.fixedWeight = fixedWeight;
M.nbDemos = nbDemos;
M.L = L;
M.Mu_d = Mu_d;
M.R_Sigma_d = R_Sigma_d;
M.G = G;
M.Mu_g = Mu_g;
M.R_Sigma_g = R_Sigma_g;
M.Mu_x = Mu_x;
M.R_Sigma_x = R_Sigma_x;
M.Demos = Demos;

metaSolver =  'matlab'; % 'pso' 'matlab' 'cmaes' 'use_existing';
nVars = 3; % number of varianbles/weights for meta optimization

switch metaSolver
    case 'cmaes'
        %% CMA-ES
        opts.LBounds = 0; opts.UBounds = 1;
        % opts.Restarts = 3;  % doubles the popsize for each restart
        doSoftConstraint = 1;
        [x, F_cmaes, E_cmaes, STOP, OUT] = cmaes('objfcn', rand(nVars,1), 1/6, opts, M, doSoftConstraint);
        plotcmaesdat
        
    case 'pso'
        %% PSO
        lb = 0*ones(1,nVars);
        ub = 1*ones(1,nVars);
        
        options = optimoptions('particleswarm','SwarmSize',2*nVars, 'Display', 'iter');
        doSoftConstraint = 1;
        
        fh = @(x)objfcn(x, M, doSoftConstraint);
        [x, fval, exitflag] = particleswarm(fh, nVars, lb, ub, options);
    case 'use_existing'
        x = [0.5 0.4910 0.0090]; % for G skill (5)
        
    case 'matlab'
        lb = 0*ones(1,nVars);
        ub = 1*ones(1,nVars);
        doSoftConstraint = 0; % no need for soft constraints since it is enforced as hard linear constraint
        
        options = optimoptions('fmincon', 'Algorithm','sqp','MaxIterations',1000); 
        
        fh = @(x)objfcn(x, M, doSoftConstraint);
        [x, fval, exitflag] = fmincon(fh, rand(nVars,1), [], [], ones(1,nVars), 1, lb, ub, [], options);
end

% output of this section is the weight between the position and shape costs

%% check the result of the meta-optimzation
w = x;     % weight
P_ = zeros( nbDims, nbNodes);
P_(1,1) = fixedWeight;
P_(2,end) = fixedWeight;

figure;
whichDemos = [1 2 3];
Sols = cell(1,length(whichDemos));
for ni = 1:length(whichDemos)
    % define the constraint
    posConstraints = [(Demos{whichDemos(ni)}(:,1)+0*rand(2,1)).' ; (Demos{whichDemos(ni)}(:,end)+0*rand(2,1)).']*fixedWeight;
    
    % CVX
    cvx_begin
    variable sol_x(nbNodes);
    variable sol_y(nbNodes);
    minimize(w(1) .*  ((R_Sigma_d * reshape((L*[sol_x sol_y] - Mu_d.').', numel(Mu_d),1)).' * (R_Sigma_d * reshape((L*[sol_x sol_y] - Mu_d.').', numel(Mu_d),1))) + ...
        w(2) .* ((R_Sigma_g * reshape((G*[sol_x sol_y] - Mu_g.').', numel(Mu_g),1)).' * (R_Sigma_g * reshape((G*[sol_x sol_y] - Mu_g.').', numel(Mu_g),1))) + ...
        w(3) .* ((R_Sigma_x * reshape(([sol_x sol_y] - Mu_x.').', numel(Mu_x),1)).' * (R_Sigma_x * reshape(([sol_x, sol_y] - Mu_x.').', numel(Mu_x),1))))
    % minimize(f([sol_x, sol_y]));
    subject to
    P_*[sol_x, sol_y] == posConstraints;
    cvx_end
    
    sol = [sol_x, sol_y];
    Sols{1,ni} = sol;
    
    % plot
    subplot(1,length(whichDemos),ni);hold on
    title('GMM- \delta');
    for ii=1:nbDemos
        plot(Demos{ii}(1,:),Demos{ii}(2,:),'color',[0.5 0.5 0.5]);
    end
    plot(sol(:,1),sol(:,2),'linewidth',2)
    plot(Mu_x(1,:),Mu_x(2,:),'--r','linewidth',2)
    bound_x = abs(max(sol_x) - min(sol_x))*0.1;
    bound_y = abs(max(sol_y) - min(sol_y))*0.1;
    axis([min(sol(:,1))-bound_x max(sol(:,1))+bound_x min(sol(:,2))-bound_y max(sol(:,2))+bound_y]);
    xticklabels([]);
    yticklabels([]);
    box on; grid on;
    ylabel('x_2','fontname','Times','fontsize',14);
    xlabel('x_1','fontname','Times','fontsize',14);
end

for ii=1:length(whichDemos); subplot(1,length(whichDemos),ii);axis auto;end
return
% save the important variables for plotting later
filenamesaved = ['skill_' num2str(kk) '_trained.mat'];
save(filenamesaved,'Demos','Gmms','Sols','w');