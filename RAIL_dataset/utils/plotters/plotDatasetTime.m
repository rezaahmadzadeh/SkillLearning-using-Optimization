function plotDatasetTime(dataset, fig)

if nargin ==1
   fig = figure;
end

plotDatasetField(dataset, 'pos', fig);

if isfield(dataset, 'vel')
    fig = figure;
    plotDatasetField(dataset, 'vel', fig);
end

if isfield(dataset, 'acc')
    fig = figure;
    plotDatasetField(dataset, 'acc', fig);
end

end

function plotDatasetField(dataset, fieldname, figureHandle)

numVars = size(eval(['dataset(1).', fieldname]),2);  

f = figure(figureHandle); hold on;
for i=1:numVars
    h = subplot(numVars,1,i);
    hold on; ylabel(['x_', num2str(i)]); xlabel('time');
    axis tight; grid on;
    
    legendVec = cell(size(dataset));
    for ind = 1:size(dataset,1)  
        pos = eval(['dataset(ind).', fieldname]);
        plot(dataset(ind).time, pos(:,i), '-', 'linewidth', 2)
        legendVec{ind} = sprintf('Demo %i', ind);
    end
    legend(legendVec);
    hold off;
end

end
