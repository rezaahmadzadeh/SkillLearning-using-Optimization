function dataset_warped = alignDataset(dataset, ind_source)

if nargin ==1
   ind_source = 1;
end

if ind_source == Inf
    dataset_warped = dataset;
    return
end

numDim = size(dataset{1}.pos,1); % number of dimensions

traj_source = [dataset{ind_source}.time', dataset{ind_source}.pos(1:numDim,:)'];

dataset_warped = cell(size(dataset));

for ii = 1:length(dataset)
    
    traj_target = [dataset{ii}.time', dataset{ii}.pos(1:numDim,:)'];
    traj_warped = DTW(traj_source', traj_target').';
    
    %         [traj_source, traj_warped, ~] = DTW_Calinon(traj_source', traj_target');
    %         traj_warped = traj_warped';
    %         traj_source = traj_source';
    
    Demo = struct();
    Demo.time = traj_warped(:,1);
    Demo.pos = traj_warped(:, 2:end);
    % dataset_warped = [dataset_warped; Demo];  % <-- Why?
    dataset_warped{ii}.time = Demo.time';
    dataset_warped{ii}.pos = Demo.pos';
%     dataset_warped(ii,1).obj = dataset{ii}.obj;
    
end

end

function [data_warped] = DTW(data_source, data_target)

nbVar = size(data_source,1);
time = data_source(1,:);


SM = simmx(data_source(2:end,:), data_target(2:end,:));
%     subplot(121)
%     imagesc(SM)
%     colormap(1-gray)

[p,q,C] = dp(SM);
%[p,q,C] = dp(1-SM);
%     hold on; plot(q,p,'r'); hold off
%     subplot(122)
%     imagesc(C)
%     hold on; plot(q,p,'r'); hold off
%

matched_ind = zeros(1, size(data_source,2));

for i = 1:length(matched_ind)
    matched_ind(i) = q(min(find(p >= i)));
end

data_warped(1,:) = time;
data_warped(2:nbVar,:) = data_target(2:end,matched_ind);

end

function M = simmx(A,B)
% M = simmx(A,B)
%    calculate a sim matrix between specgram-like feature matrices A and B.
%    size(M) = [size(A,2) size(B,2)]; A and B have same #rows.
% 2003-03-15 dpwe@ee.columbia.edu

% Copyright (c) 2003 Dan Ellis <dpwe@ee.columbia.edu>
% released under GPL - see file COPYRIGHT

% EA = sqrt(sum(A.^2));
% EB = sqrt(sum(B.^2));

%ncA = size(A,2);
%ncB = size(B,2);
%M = zeros(ncA, ncB);
%for i = 1:ncA
%  for j = 1:ncB
%    % normalized inner product i.e. cos(angle between vectors)
%    M(i,j) = (A(:,i)'*B(:,j))/(EA(i)*EB(j));
%  end
%end

% this is 10x faster
M = pdist2(A',B');
%M = (A'*B)./(EA'*EB);
end

function [p,q,D] = dp(M)
% [p,q] = dp(M) 
%    Use dynamic programming to find a min-cost path through matrix M.
%    Return state sequence in p,q
% 2003-03-15 dpwe@ee.columbia.edu

% Copyright (c) 2003 Dan Ellis <dpwe@ee.columbia.edu>
% released under GPL - see file COPYRIGHT

[r,c] = size(M);

% costs
D = zeros(r+1, c+1);
D(1,:) = NaN;
D(:,1) = NaN;
D(1,1) = 0;
D(2:(r+1), 2:(c+1)) = M;

% traceback
phi = zeros(r,c);

for i = 1:r; 
  for j = 1:c;
    [dmax, tb] = min([D(i, j), D(i, j+1), D(i+1, j)]);
    D(i+1,j+1) = D(i+1,j+1)+dmax;
    phi(i,j) = tb;
  end
end

% Traceback from top left
i = r; 
j = c;
p = i;
q = j;
while i > 1 & j > 1
  tb = phi(i,j);
  if (tb == 1)
    i = i-1;
    j = j-1;
  elseif (tb == 2)
    i = i-1;
  elseif (tb == 3)
    j = j-1;
  else    
    error;
  end
  p = [i,p];
  q = [j,q];
end

% Strip off the edges of the D matrix before returning
 D = D(2:(r+1),2:(c+1));
end
