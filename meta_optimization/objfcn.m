function J = objfcn(w, M, doSoftConstraint)
% w is our lambda, weight between the costs
% M is a structure containing all the required arguments

% unpacking M
nbDims = M.nbDims;
nbNodes = M.nbNodes;
fixedWeight = M.fixedWeight;
nbDemos = M.nbDemos;
L = M.L;
Mu_d = M.Mu_d;
Mu_x = M.Mu_x;
R_Sigma_d = M.R_Sigma_d;
R_Sigma_x = M.R_Sigma_x;
Demos = M.Demos;

P_ = zeros( nbDims, nbNodes);
P_(1,1) = fixedWeight;
P_(2,end) = fixedWeight;
% if w > 1, w=1; end
% if w < 0, w=0; end
disp(['Weights: ' num2str(w(1)) ' , ' num2str(w(2))]);
Sols = cell(1,nbDemos);
% get the initial point form another demo

for ni = 1:nbDemos
    % define the constraint
    G = [(Demos{ni}(:,1)+0*rand(2,1)).' ; (Demos{ni}(:,end)+0*rand(2,1)).']*fixedWeight;
    
    % CVX
    cvx_begin quiet
    variable sol_x(nbNodes);
    variable sol_y(nbNodes);
    minimize(w(1) .*  ((R_Sigma_d * reshape((L*[sol_x sol_y] - Mu_d.').', numel(Mu_d),1)).' * (R_Sigma_d * reshape((L*[sol_x sol_y] - Mu_d.').', numel(Mu_d),1))) + ...
        w(2) .* ((R_Sigma_x * reshape(([sol_x sol_y] - Mu_x.').', numel(Mu_x),1)).' * (R_Sigma_x * reshape(([sol_x, sol_y] - Mu_x.').', numel(Mu_x),1))))
    % minimize(f([sol_x, sol_y]));
    subject to
    P_*[sol_x, sol_y] == G;
    cvx_end
    
    sol = [sol_x, sol_y];
    Sols{1,ni} = sol;
end

J = 0;
for ii=1:size(Demos,2)
    J = J + (sum(sum((Sols{ii} - Demos{ii}.').^2)));
end

if doSoftConstraint
    J = J + 1e10*abs((1-(w(1)+w(2)))); % soft constraint to normalize w1 and w2
end
end
