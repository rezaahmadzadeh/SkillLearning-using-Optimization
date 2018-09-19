function datasetChangePointed = findCurvatureChangePointsDataset(dataset, numChangepoints, endCut)

if nargin==2
    endCut = 100;
end

datasetChangePointed = dataset;

for ii=1:length(dataset)
%     dr = dataset(ii).vel;
%     ddr = dataset(ii).acc;
    nPts = size(dataset(ii).pos,1);
    [~, dr, ddr] = diffArcLength(dataset(ii).pos);
    
    K = findCurvatureTrajectory(dr, ddr);
    
    [psor,lsor] = findpeaks(K, 'SortStr','descend');
    
    lsor(lsor<endCut | lsor>(nPts-endCut)) = []; 
    
    if length(lsor)>=numChangepoints
        lsor = lsor(1:numChangepoints);
        lsor = sort(lsor);
        segs = [1, lsor, nPts];
    else
        lsor = [];
        segs = [];
    end
    
    datasetChangePointed(ii).changeIdx = lsor;
    datasetChangePointed(ii).segIdx = segs;
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function K = findCurvatureTrajectory(dr, ddr)

for n=1:size(dr,1)
   
    cross_prod = norm(cross(dr(n,:)', ddr(n,:)'));
    
    K(n) = cross_prod/(norm(dr(n,:)'))^3;
    
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [arcLength, dTraj, d2Traj] = diffArcLength(traj)

numDim = size(traj,2);
numPts = size(traj,1);


dTraj = zeros(numPts, numDim);
d2Traj = zeros(numPts, numDim);

% compute the chord linear arclengths
seglen = sqrt(sum(diff(traj,[],1).^2,2));
chordlen = seglen;

% compute the splines
spl = cell(1,numDim);
spld = spl;
spld2 = spl;

diffarray = [3 0 0;0 2 0;0 0 1;0 0 0];
diff2array = [2 0;0 1;0 0];

for i = 1:numDim
    spl{i} = spline([0;cumsum(chordlen)],traj(:,i));
    nc = numel(spl{i}.coefs);
    
    if nc < 4
        % just pretend it has cubic segments
        spl{i}.coefs = [zeros(1,4-nc),spl{i}.coefs];
        spl{i}.order = 4;
    end

    % and now differentiate them
    xp = spl{i};
    xp.coefs = xp.coefs*diffarray;
    xp.order = 3;
    spld{i} = xp;
   
end

% numerical integration along the curve
polyarray = zeros(numDim,3);

for i = 1:spl{1}.pieces
  % extract polynomials for the derivatives
  for j = 1:numDim
    polyarray(j,:) = spld{j}.coefs(i,:);
  end
  
  % integrate the arclength for the i'th segment
  % using quadgk for the integral. I could have
  % done this part with an ode solver too.
  seglen(i) = quadgk(@(t) segkernel(t),0,chordlen(i));
end

% and sum the segments
arcLength = [0;cumsum(seglen)];

for i = 1:numDim
    spl{i} = spline(arcLength,traj(:,i));
    nc = numel(spl{i}.coefs);
    
    if nc < 4
        % just pretend it has cubic segments
        spl{i}.coefs = [zeros(1,4-nc),spl{i}.coefs];
        spl{i}.order = 4;
    end

    % and now differentiate them
    xp = spl{i};
    xp.coefs = xp.coefs*diffarray;
    xp.order = 3;
    spld{i} = xp;
    
    xp = spld{i};
    xp.coefs = xp.coefs*diff2array;
    xp.order = 2;
    spld2{i} = xp;
end

% numerical integration along the curve
polyarray = zeros(numDim,3);
polyarray2 = zeros(numDim,2);

for i = 1:spl{1}.pieces
  % extract polynomials for the derivatives
  for j = 1:numDim
    polyarray(j,:) = spld{j}.coefs(i,:);
    dTraj(i+1,j) = polyval(polyarray(j,:),seglen(i));
    polyarray2(j,:) = spld2{j}.coefs(i,:);
    d2Traj(i+1,j) = polyval(polyarray2(j,:),seglen(i));
  end
end


for j=1:numDim
    polyarray(j,:) = spld{j}.coefs(1,:);
    polyarray2(j,:) = spld2{j}.coefs(1,:);
    dTraj(1,j) = polyval(polyarray(j,:),0);
    d2Traj(1,j) = polyval(polyarray2(j,:),0);
end

% ==========================
%   end main function
% ==========================
%   begin nested functions
% ==========================
  function val = segkernel(t)
    % sqrt((dx/dt)^2 + (dy/dt)^2)
    
    val = zeros(size(t));
    for k = 1:numDim
      val = val + polyval(polyarray(k,:),t).^2;
    end
    val = sqrt(val);
    
  end % function segkernel

end % function arclength