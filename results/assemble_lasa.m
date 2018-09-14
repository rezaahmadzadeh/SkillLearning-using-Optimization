clc,clear,%close all

c(1)= openfig('angle_trained.fig','invisible');
c(2)=openfig('BendedLine_trained.fig','invisible');
c(3) = openfig('CShape_trained.fig','invisible');
c(4) = openfig('DoubleBendedLine_trained.fig','invisible');
c(5) = openfig('GShape_trained.fig','invisible');
c(6) = openfig('heee_trained.fig','invisible');
c(7) = openfig('JShape_trained.fig','invisible');
c(8) = openfig('Khamesh_trained.fig','invisible');
c(9) = openfig('Leaf_1_trained.fig','invisible');
c(10) = openfig('Leaf_2_trained.fig','invisible');
c(11)= openfig('Line_trained.fig','invisible');
c(12) = openfig('LShape_trained.fig','invisible');
c(13) = openfig('NShape_trained.fig','invisible');
c(14) = openfig('PShape_trained.fig','invisible');
c(15) = openfig('RShape_trained.fig','invisible');
c(16) = openfig('Saeghe_trained.fig','invisible');
c(17) = openfig('Sharpc_trained.fig','invisible');
c(18) = openfig('Sine_trained.fig','invisible');
c(19) = openfig('Snake_trained.fig','invisible');
c(20) = openfig('Spoon_trained.fig','invisible');
c(21) = openfig('Sshape_trained.fig','invisible');
c(22) = openfig('Trapezoid_trained.fig','invisible');
c(23) = openfig('Worm_trained.fig','invisible');
c(24) = openfig('WShape_trained.fig','invisible');
c(25) = openfig('Zshape_trained.fig','invisible');


for ii=1:25
set(findobj(c(ii), 'Type', 'Line', 'color', 'b'), 'LineWidth', 1);
set(findobj(c(ii), 'Type', 'Line', 'color', 'r'), 'LineWidth', 2);
% set(findobj(c1, 'Type', 'Line', 'color', 'r'), 'Linestyle', '-');
end

% return
% H = figure;
H = figure('units','normalized','outerposition',[0 0 1 1]);

for ii=1:25
    h(ii) = subplot(5,5,ii);
    copyobj(allchild(get(c(ii),'CurrentAxes')),h(ii));
    set(h(ii),'xtick',[]);set(h(ii),'ytick',[]);
    subplot(h(ii));box on; axis square
end

n = 0.1;
for ii=1:25
    subplot(5,5,ii);
    X = xlim;
    s = n*abs((X(2) - X(1)));
    X(1) = X(1) - s;
    X(2) = X(2) + s;
    xlim([X(1) X(2)]);
    
    Y = ylim;
    s = n*abs((Y(2) - Y(1)));
    Y(1) = Y(1) - s;
    Y(2) = Y(2) + s;
    ylim([Y(1) Y(2)]);    
end


hd = 0.1;
vd = 0.04;
for ii=[2,7,12,17,22]
p = get(h(ii), 'position');
p(1) = p(1) - hd;
set(h(ii),'position',p);
end

for ii=[3,8,13,18,23]
p = get(h(ii), 'position');
p(1) = p(1) - 2*hd;
set(h(ii),'position',p);
end

for ii=[4,9,14,19,24]
p = get(h(ii), 'position');
p(1) = p(1) - 3*hd;
set(h(ii),'position',p);
end

for ii=[5,10,15,20,25]
p = get(h(ii), 'position');
p(1) = p(1) - 4*hd;
set(h(ii),'position',p);
end

for ii = 6:10
p = get(h(ii), 'position');
p(2) = p(2) + vd;
set(h(ii),'position',p);
end


for ii = 11:15
p = get(h(ii), 'position');
p(2) = p(2) + 2*vd;
set(h(ii),'position',p);
end

for ii = 16:20
p = get(h(ii), 'position');
p(2) = p(2) + 3*vd;
set(h(ii),'position',p);
end


for ii = 21:25
p = get(h(ii), 'position');
p(2) = p(2) + 4*vd;
set(h(ii),'position',p);
end

    

