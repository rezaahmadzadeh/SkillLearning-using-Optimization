function [evaluationResults, repros] = compare_with_baslines(resultsFile,doPlot)
%COMPARE_WITH_BASLINES compares automated cost balancing algorithm with
%
%   four baselines: (1) position-coordinate only; (2) tangent
%   coordinate only; (3) laplacian-coordinate only; (4) uniform weighting
%   of all three coordinates
%   

%% Load results
dataObj = load(resultsFile);
M = dataObj.M;
demos = dataObj.Demos;

%% Unpack M

nbDims = M.nbDims;
nbNodes = M.nbNodes;
fixedWeight = M.fixedWeight;
nbDemos = M.nbDemos;
L = M.L;
G = M.G;
Mu_d = M.Mu_d;
Mu_g = M.Mu_g;
Mu_x = M.Mu_x;
R_Sigma_d = M.R_Sigma_d;
R_Sigma_g = M.R_Sigma_g;
R_Sigma_x = M.R_Sigma_x;
Demos = M.Demos;
viaPoints = M.viaPoints;
viaPointsTime = M.viaPointsTime;
doConstraintIntialPoint = M.doConstraintIntialPoint;
doConstraintEndPoint = M.doConstraintEndPoint;

%% Define P_

nbViaPoints = length(viaPointsTime);
nbConstraintPoints = nbViaPoints + doConstraintIntialPoint + doConstraintEndPoint;

P_ = zeros((nbConstraintPoints), nbNodes);

P_index = 1;
if(doConstraintIntialPoint)
    P_(P_index,1) = fixedWeight; % initial point
    P_index = P_index + 1;
end

if(doConstraintEndPoint)
    P_(P_index,end) = fixedWeight; % end point
    P_index = P_index + 1;
end

for i = 1:nbViaPoints
    P_(P_index,viaPointsTime(i)) = fixedWeight;
    
    P_index = P_index + 1;
end

%% Generate reproductions

repros{1} = cell(1,nbDemos);
repros{2} = cell(1,nbDemos);
repros{3} = cell(1,nbDemos);
repros{4} = cell(1,nbDemos);
repros{5} = cell(1,nbDemos);

for ni = 1:nbDemos
    % define the constraint
    posConstraints = [(Demos{ni}(:,1)+0*rand(nbDims,1)).' ; (Demos{ni}(:,end)+0*rand(nbDims,1)).']*fixedWeight;
    
    if ~isempty(viaPoints) % hard coded via points enforced independent of demonstrations
        posConstraints = [posConstraints; viaPoints.'.*fixedWeight];
    else
        if ~isempty(viaPointsTime)
            for i = 1:nbViaPoints
                posConstraints = [posConstraints; Demos{ni}(:,viaPointsTime(i)).'*fixedWeight];
            end
        end
    end
    
    % position-only coordinate
    w = [0 0 1];
    scalingFactors = ones(size(w));
    
    if nbDims == 2
        cvx_begin quiet
        variable sol_x(nbNodes);
        variable sol_y(nbNodes);
        minimize(w(1) .*  ((R_Sigma_d * reshape((L*[sol_x sol_y] - Mu_d.').', numel(Mu_d),1)).' * (R_Sigma_d * reshape((L*[sol_x sol_y] - Mu_d.').', numel(Mu_d),1)))./scalingFactors(1) + ...
            w(2) .* ((R_Sigma_g * reshape((G*[sol_x sol_y] - Mu_g.').', numel(Mu_g),1)).' * (R_Sigma_g * reshape((G*[sol_x sol_y] - Mu_g.').', numel(Mu_g),1)))./scalingFactors(2) + ...
            w(3) .* ((R_Sigma_x * reshape(([sol_x sol_y] - Mu_x.').', numel(Mu_x),1)).' * (R_Sigma_x * reshape(([sol_x, sol_y] - Mu_x.').', numel(Mu_x),1)))./scalingFactors(3))
            % minimize(f([sol_x, sol_y]));
        subject to
        P_*[sol_x, sol_y] == posConstraints;
        cvx_end
        sol = [sol_x, sol_y];
    else
        if nbDims == 3
            cvx_begin quiet
            variable sol_x(nbNodes);
            variable sol_y(nbNodes);
            variable sol_z(nbNodes);
            minimize(w(1) .*  ((R_Sigma_d * reshape((L*[sol_x sol_y sol_z] - Mu_d.').', numel(Mu_d),1)).' * (R_Sigma_d * reshape((L*[sol_x sol_y sol_z] - Mu_d.').', numel(Mu_d),1)))./scalingFactors(1) + ...
                w(2) .* ((R_Sigma_g * reshape((G*[sol_x sol_y sol_z] - Mu_g.').', numel(Mu_g),1)).' * (R_Sigma_g * reshape((G*[sol_x sol_y sol_z] - Mu_g.').', numel(Mu_g),1)))./scalingFactors(2) + ...
                w(3) .* ((R_Sigma_x * reshape(([sol_x sol_y sol_z] - Mu_x.').', numel(Mu_x),1)).' * (R_Sigma_x * reshape(([sol_x, sol_y sol_z] - Mu_x.').', numel(Mu_x),1)))./scalingFactors(3))            
            subject to
            P_*[sol_x, sol_y sol_z] == posConstraints;
            cvx_end
            sol = [sol_x, sol_y sol_z];
        else
            error("The current version of the software can only handle 2 and 3 dimensional spaces!")
        end
    end
    
    repros{1}{1,ni} = sol.';
    clear sol;
    
    % tangent-only coordinate
    w = [0 1 0];
    scalingFactors = ones(size(w));
    
    if nbDims == 2
        cvx_begin quiet
        variable sol_x(nbNodes);
        variable sol_y(nbNodes);
        minimize(w(1) .*  ((R_Sigma_d * reshape((L*[sol_x sol_y] - Mu_d.').', numel(Mu_d),1)).' * (R_Sigma_d * reshape((L*[sol_x sol_y] - Mu_d.').', numel(Mu_d),1)))./scalingFactors(1) + ...
            w(2) .* ((R_Sigma_g * reshape((G*[sol_x sol_y] - Mu_g.').', numel(Mu_g),1)).' * (R_Sigma_g * reshape((G*[sol_x sol_y] - Mu_g.').', numel(Mu_g),1)))./scalingFactors(2) + ...
            w(3) .* ((R_Sigma_x * reshape(([sol_x sol_y] - Mu_x.').', numel(Mu_x),1)).' * (R_Sigma_x * reshape(([sol_x, sol_y] - Mu_x.').', numel(Mu_x),1)))./scalingFactors(3))
            % minimize(f([sol_x, sol_y]));
        subject to
        P_*[sol_x, sol_y] == posConstraints;
        cvx_end
        sol = [sol_x, sol_y];
    else
        if nbDims == 3
            cvx_begin quiet
            variable sol_x(nbNodes);
            variable sol_y(nbNodes);
            variable sol_z(nbNodes);
            minimize(w(1) .*  ((R_Sigma_d * reshape((L*[sol_x sol_y sol_z] - Mu_d.').', numel(Mu_d),1)).' * (R_Sigma_d * reshape((L*[sol_x sol_y sol_z] - Mu_d.').', numel(Mu_d),1)))./scalingFactors(1) + ...
                w(2) .* ((R_Sigma_g * reshape((G*[sol_x sol_y sol_z] - Mu_g.').', numel(Mu_g),1)).' * (R_Sigma_g * reshape((G*[sol_x sol_y sol_z] - Mu_g.').', numel(Mu_g),1)))./scalingFactors(2) + ...
                w(3) .* ((R_Sigma_x * reshape(([sol_x sol_y sol_z] - Mu_x.').', numel(Mu_x),1)).' * (R_Sigma_x * reshape(([sol_x, sol_y sol_z] - Mu_x.').', numel(Mu_x),1)))./scalingFactors(3))            
            subject to
            P_*[sol_x, sol_y sol_z] == posConstraints;
            cvx_end
            sol = [sol_x, sol_y sol_z];
        else
            error("The current version of the software can only handle 2 and 3 dimensional spaces!")
        end
    end
    
    repros{2}{1,ni} = sol.';
    clear sol;
    
    % laplace-only coordinate
    w = [1 0 0];
    scalingFactors = ones(size(w));
    
    if nbDims == 2
        cvx_begin quiet
        variable sol_x(nbNodes);
        variable sol_y(nbNodes);
        minimize(w(1) .*  ((R_Sigma_d * reshape((L*[sol_x sol_y] - Mu_d.').', numel(Mu_d),1)).' * (R_Sigma_d * reshape((L*[sol_x sol_y] - Mu_d.').', numel(Mu_d),1)))./scalingFactors(1) + ...
            w(2) .* ((R_Sigma_g * reshape((G*[sol_x sol_y] - Mu_g.').', numel(Mu_g),1)).' * (R_Sigma_g * reshape((G*[sol_x sol_y] - Mu_g.').', numel(Mu_g),1)))./scalingFactors(2) + ...
            w(3) .* ((R_Sigma_x * reshape(([sol_x sol_y] - Mu_x.').', numel(Mu_x),1)).' * (R_Sigma_x * reshape(([sol_x, sol_y] - Mu_x.').', numel(Mu_x),1)))./scalingFactors(3))
            % minimize(f([sol_x, sol_y]));
        subject to
        P_*[sol_x, sol_y] == posConstraints;
        cvx_end
        sol = [sol_x, sol_y];
    else
        if nbDims == 3
            cvx_begin quiet
            variable sol_x(nbNodes);
            variable sol_y(nbNodes);
            variable sol_z(nbNodes);
            minimize(w(1) .*  ((R_Sigma_d * reshape((L*[sol_x sol_y sol_z] - Mu_d.').', numel(Mu_d),1)).' * (R_Sigma_d * reshape((L*[sol_x sol_y sol_z] - Mu_d.').', numel(Mu_d),1)))./scalingFactors(1) + ...
                w(2) .* ((R_Sigma_g * reshape((G*[sol_x sol_y sol_z] - Mu_g.').', numel(Mu_g),1)).' * (R_Sigma_g * reshape((G*[sol_x sol_y sol_z] - Mu_g.').', numel(Mu_g),1)))./scalingFactors(2) + ...
                w(3) .* ((R_Sigma_x * reshape(([sol_x sol_y sol_z] - Mu_x.').', numel(Mu_x),1)).' * (R_Sigma_x * reshape(([sol_x, sol_y sol_z] - Mu_x.').', numel(Mu_x),1)))./scalingFactors(3))            
            subject to
            P_*[sol_x, sol_y sol_z] == posConstraints;
            cvx_end
            sol = [sol_x, sol_y sol_z];
        else
            error("The current version of the software can only handle 2 and 3 dimensional spaces!")
        end
    end
    
    repros{3}{1,ni} = sol.';
    clear sol;
    
    % uniform weighting
    w = [1/3 1/3 1/3];
    scalingFactors = ones(size(w));
    
    if nbDims == 2
        cvx_begin quiet
        variable sol_x(nbNodes);
        variable sol_y(nbNodes);
        minimize(w(1) .*  ((R_Sigma_d * reshape((L*[sol_x sol_y] - Mu_d.').', numel(Mu_d),1)).' * (R_Sigma_d * reshape((L*[sol_x sol_y] - Mu_d.').', numel(Mu_d),1)))./scalingFactors(1) + ...
            w(2) .* ((R_Sigma_g * reshape((G*[sol_x sol_y] - Mu_g.').', numel(Mu_g),1)).' * (R_Sigma_g * reshape((G*[sol_x sol_y] - Mu_g.').', numel(Mu_g),1)))./scalingFactors(2) + ...
            w(3) .* ((R_Sigma_x * reshape(([sol_x sol_y] - Mu_x.').', numel(Mu_x),1)).' * (R_Sigma_x * reshape(([sol_x, sol_y] - Mu_x.').', numel(Mu_x),1)))./scalingFactors(3))
            % minimize(f([sol_x, sol_y]));
        subject to
        P_*[sol_x, sol_y] == posConstraints;
        cvx_end
        sol = [sol_x, sol_y];
    else
        if nbDims == 3
            cvx_begin quiet
            variable sol_x(nbNodes);
            variable sol_y(nbNodes);
            variable sol_z(nbNodes);
            minimize(w(1) .*  ((R_Sigma_d * reshape((L*[sol_x sol_y sol_z] - Mu_d.').', numel(Mu_d),1)).' * (R_Sigma_d * reshape((L*[sol_x sol_y sol_z] - Mu_d.').', numel(Mu_d),1)))./scalingFactors(1) + ...
                w(2) .* ((R_Sigma_g * reshape((G*[sol_x sol_y sol_z] - Mu_g.').', numel(Mu_g),1)).' * (R_Sigma_g * reshape((G*[sol_x sol_y sol_z] - Mu_g.').', numel(Mu_g),1)))./scalingFactors(2) + ...
                w(3) .* ((R_Sigma_x * reshape(([sol_x sol_y sol_z] - Mu_x.').', numel(Mu_x),1)).' * (R_Sigma_x * reshape(([sol_x, sol_y sol_z] - Mu_x.').', numel(Mu_x),1)))./scalingFactors(3))            
            subject to
            P_*[sol_x, sol_y sol_z] == posConstraints;
            cvx_end
            sol = [sol_x, sol_y sol_z];
        else
            error("The current version of the software can only handle 2 and 3 dimensional spaces!")
        end
    end
    
    repros{4}{1,ni} = sol.';
    clear sol;
end

%% optimal weighting (ours)

% w = dataObj.w;
% scalingFactors = M.scalingFactors; 

repros{5} = dataObj.Sols;

%% evaluate the reproductions

evaluationResults{1}.algoName = 'position-only';
evaluationResults{2}.algoName = 'tangent-only';
evaluationResults{3}.algoName = 'laplace-only';
evaluationResults{4}.algoName = 'uniform weighting';
evaluationResults{5}.algoName = 'optimal weighting';

for i = 1:5
    evaluationResults{i}.performanceMeasures = evaluate_reproductions(demos,repros{i});
    if doPlot
        figure(1)
        title('Swept Error Area') 
        plot(evaluationResults{i}.performanceMeasures.SEA.list,'lineWidth',2)
        legend(evaluationResults{1}.algoName,evaluationResults{2}.algoName,evaluationResults{3}.algoName,evaluationResults{4}.algoName,evaluationResults{5}.algoName)
        hold on
        
        figure(2)
        title('Sum of Squared Errors')
        plot(evaluationResults{i}.performanceMeasures.SSE.list,'lineWidth',2)
        legend(evaluationResults{1}.algoName,evaluationResults{2}.algoName,evaluationResults{3}.algoName,evaluationResults{4}.algoName,evaluationResults{5}.algoName)
        hold on
        
        figure(3)
        title('Dynamic Time Warping Distance')
        plot(evaluationResults{i}.performanceMeasures.DTWD.list,'LineWidth',2)
        legend(evaluationResults{1}.algoName,evaluationResults{2}.algoName,evaluationResults{3}.algoName,evaluationResults{4}.algoName,evaluationResults{5}.algoName)
        hold on
        
        figure(4)
        title('Hausdorff Distance')
        plot(evaluationResults{i}.performanceMeasures.HD.list,'LineWidth',2)
        legend(evaluationResults{1}.algoName,evaluationResults{2}.algoName,evaluationResults{3}.algoName,evaluationResults{4}.algoName,evaluationResults{5}.algoName)
        hold on
        
        figure(5)
        title('Frechet Distance')
        plot(evaluationResults{i}.performanceMeasures.FD.list,'LineWidth',2)
        legend(evaluationResults{1}.algoName,evaluationResults{2}.algoName,evaluationResults{3}.algoName,evaluationResults{4}.algoName,evaluationResults{5}.algoName)
        hold on
    end
end

