function dataset_diffed = differentiateDataset(dataset, method, idx)
% method can either be spline or diff or pchip
% NOTE: Finite differencing starts and end at zero vel/acc which is
% desirable!!

if nargin == 1
    method = 'diff';
    idx = [1,2,3];
elseif nargin == 2
    idx = [1,2,3];
end

dataset_diffed = dataset;

if strcmp(method, 'diff')
    dt = dataset(1).time(2) - dataset(1).time(1);
    
    for j = 1:size(dataset,1)
        time = dataset(j).time;
        %time = [0; [dt + time]];
        
        pos = dataset(j).pos(:,idx);
        %pos = [repmat(pos(1,:),1,1); pos];
        
        numVars = size(pos,2);
        numPts = size(time,1);
        
        pos_diffed = zeros(numPts, numVars);
        vel_diffed = zeros(numPts, numVars);
        acc_diffed = zeros(numPts, numVars);
        
        for i = 1:numVars
            pos_diffed(:,i) = pos(:,i);
            vel_diffed(1:(numPts-1),i) = diff(pos(:,i))/dt;
            acc_diffed(1:(numPts-1),i) = diff(vel_diffed(:,i))/dt;
            
        end
        dataset_diffed(j).time = time;
        dataset_diffed(j).pos = pos_diffed;
        dataset_diffed(j).vel = vel_diffed;
        dataset_diffed(j).acc = acc_diffed;

    end
% For pchip or spline    
else
    for j=1:size(dataset,1)
        time = dataset(j).time;
        pos = dataset(j).pos;
        
        numVars = size(pos,2);
        numPts = size(time,1);
        
        pos_diffed = zeros(numPts, numVars);
        vel_diffed = zeros(numPts, numVars);
        acc_diffed = zeros(numPts, numVars);
        
        for i=1:numVars
            %spl = eval([method, '(time, pos(:,i))']);
            %spl = spaps(time, pos(:,i), 1e-4);
            spl = spline(time, [0; pos(:,i); 0]);
            pos_diffed(:,i) = fnval(spl, time);
            vel_diffed(:,i) = fnval(fnder(spl,1),time);
            acc_diffed(:,i) = fnval(fnder(spl,2),time);
        end
        
        dataset_diffed(j).time = time;
        dataset_diffed(j).pos = pos_diffed;
        dataset_diffed(j).vel = vel_diffed;
        dataset_diffed(j).acc = acc_diffed;
        
    end
end




end