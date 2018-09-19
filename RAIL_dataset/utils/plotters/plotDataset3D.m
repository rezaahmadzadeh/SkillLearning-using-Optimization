function figHandle = plotDataset3D(dataset, figHandle, plotColor, legendFlag)

if nargin ==1
    figHandle = figure;
    plotColor = 'b';
    legendFlag = 1;
elseif nargin ==2
    plotColor = 'b';
    legendFlag =1;
elseif nargin == 3
    legendFlag =1;
end

set(0, 'currentfigure', figHandle);
hold on;

n = numel(dataset);
legendVec = cell(n*2,1);
%C = colorcube(size(dataset,1));

for ind = 1:n
    %     plot3(dataset(ind).pos(:,1), dataset(ind).pos(:,2), dataset(ind).pos(:,3),'linewidth', 2, 'color', C(ind,:));
    plot3(dataset(ind).pos(:,1), dataset(ind).pos(:,2), dataset(ind).pos(:,3),'linewidth', 2, 'color', plotColor);
    legendVec{ind} = sprintf('Demo %i', ind);
end

for ind = 1:n
    plot3(dataset(ind).obj(1,1), dataset(ind).obj(1,2), dataset(ind).obj(1,3), 'ro');
    legendVec{ind+n} = sprintf('Obj %i', ind);
end
if legendFlag
    legend(legendVec);
end
xlabel('x_1'); ylabel('x_2'); zlabel('x_3');
grid on;box on;
axis equal;

end

