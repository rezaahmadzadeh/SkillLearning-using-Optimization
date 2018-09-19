function [figHandles] = plotGPTime(time, Mu, Sigma, numDim, numStates, numPts, color)
% Plots a Gaussian over time
% Inputs:
%       time: a vector of time (1 x N)
%       mu: mean vector (D x N) where N is the num_pts
%       Sigma: cell of covariances (D x D x N)
%       stateInd: the index of the selected state 
%       color: color of the mean in normalized RGB vector
%       figHandle: handle of the figure to be plotted on

lightcolor = color + [0.6,0.6,0.6];
lightcolor(find(lightcolor>1.0)) = 1.0;
N = sqrt(2);

figHandles = cell(1,numStates);

for ii =1:numStates
    figHandles{ii} = figure; hold on;
    
    shiftInd = numStates*numDim*(0:(numPts-1)) + numDim*(ii - 1);
    stateInd = repelem(shiftInd,numDim) + repmat(1:numDim,1,numPts);

    stateMu = Mu(stateInd,:);
    stateSigma = Sigma(stateInd,stateInd);
    
    for jj = 1:numDim
        subplot(numDim, 1, jj); hold on; grid on; box on; axis tight
        
        dimMu = stateMu(jj:numDim:end, :);
        dimSigma = stateSigma(jj:numDim:end, jj:numDim:end);
        
        ymax = zeros(1, numPts);
        ymin = zeros(1, numPts);
        
        for n=1:numPts
             ymax(n) = dimMu(n) + sqrt(N^2.*dimSigma(n, n));
             ymin(n) = dimMu(n) - sqrt(N^2.*dimSigma(n, n));
        end
        patch([time(1:end)', time(end:-1:1)'], [ymax(1:end) ymin(end:-1:1)], lightcolor, 'LineStyle', 'none', 'FaceAlpha', 0.6);
        plot(time, dimMu, 'lineWidth', 2.5, 'color', color);
        hold off;
                
    end
end



end



