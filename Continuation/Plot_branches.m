function h = Plot_branches(parameter, a, e, E, p, varargin)
% This function allows one to plot a varying number of arrays stored in
% varargin
h = figure;
sgtitle(['Plots of the variables as a function of ', parameter, '. The initial parameter values are (a,e,E,p) = (', num2str(a), ', ', num2str(e), ', ', num2str(E), ', ', num2str(p), ').']);

vars = {'t_{0}', 'x_{0}', 'v_{0}', 'w_{0}', 't_{1}', 'x_{1}', 'v_{1}', 'w_{1}', parameter};

numPlotInputs = nargin-5; % number of arrays to plot
periodic_orbits = varargin; % n x (9xm) array where m=numPlotsInputs
% Defining the subplots
for i0=1:8
    ax(i0) = subplot(2,4,i0);        
    title(vars(i0));
    xtickformat(ax(i0), '%.2g');
    ytickformat(ax(i0), '%.2g');
    xlim auto
    ylim auto
    hold on
    for i1=1:numPlotInputs
        array = periodic_orbits{i1};
        x = array(:,9);
        y = array(:,i0);
        plot(ax(i0), x, y, 'Linewidth', 3/4, 'Color', 'black');
    end
    hold off
end
end
