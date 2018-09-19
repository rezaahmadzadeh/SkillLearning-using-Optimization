function dataset_diffed =  smoothSplineResampleDataset(dataset, num_samples, tol)

dataset_diffed = dataset;

for jj=1:numel(dataset)
    time = dataset(jj).time;
    pos = dataset(jj).pos(:,1:3);
    
    time_diffed = linspace(time(1), time(end), num_samples).';
    
    num_vars = size(pos,2);
    
    pos_diffed = zeros(num_samples, num_vars);
    vel_diffed = zeros(num_samples, num_vars);
    acc_diffed = zeros(num_samples, num_vars);
    
    for ii=1:num_vars
        spl = spaps(time, pos(:,ii), tol);
%         temp = fnval(spl, time);
%         spl = spline(time, [0; temp; 0]);
        
        
%         spl = spaps(time, pos(:,i), tol);
        pos_diffed(:,ii) = fnval(spl, time_diffed);
        vel_diffed(:,ii) = fnval(fnder(spl,1),time_diffed);
        acc_diffed(:,ii) = fnval(fnder(spl,2),time_diffed);
    end
    
    
    dataset_diffed(jj).time = time_diffed;
    dataset_diffed(jj).pos = pos_diffed;
    dataset_diffed(jj).vel = vel_diffed;
    dataset_diffed(jj).acc = acc_diffed;
    
end

end