function [repro, expSigma] = reproGMM(M)
Mu = M.Mu;
Priors = M.Priors;
Sigma = M.Sigma;
%% Use GMR to  reproduce
nbVar = size(M.Mu,1);
repro(1,:) = M.t;
[repro(2:nbVar,:), expSigma] = GMR(Priors, Mu, Sigma,  repro(1,:), 1, 2:nbVar);
end
