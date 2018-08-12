function M = encodeGMM(data, numPoints, nbStates)
M = struct();
%% Training of GMM by EM algorithm
[Priors, Mu, Sigma] = EM_init_kmeans(data, nbStates);
[Priors, Mu, Sigma] = EM(data, Priors, Mu, Sigma);

M.Priors = Priors;
M.Mu = Mu;
M.Sigma = Sigma;
M.t = linspace(min(data(1,:)), max(data(1,:)), numPoints);
end



