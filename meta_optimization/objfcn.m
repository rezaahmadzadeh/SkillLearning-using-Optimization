function J = objfcn(w, M, doSoftConstraint)
% w is our lambda, weight between the costs
% M is a structure containing all the required arguments

% unpacking M
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
scalingFactors = M.scalingFactors;
viaPoints = M.viaPoints;
viaPointsTime = M.viaPointsTime;
doConstraintIntialPoint = M.doConstraintIntialPoint;
doConstraintEndPoint = M.doConstraintEndPoint;

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

% if w > 1, w=1; end
% if w < 0, w=0; end
disp(['Weights: ' num2str(reshape(w,[1,length(w)]))]);

Sols = cell(1,nbDemos);
% get the initial point form another demo

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
    
    % CVX
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
    
    Sols{1,ni} = sol;
end

J = 0;
for ii=1:size(Demos,2)
    J = J + (sum(sum((Sols{ii} - Demos{ii}.').^2)));
end

if doSoftConstraint
    J = J + 1e10*abs((1-sum(w))); % soft constraint to normalize the weights
end
end
