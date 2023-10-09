%%%%% Compare the velocity of each crustal types (and the subclasses) at 4 depths
%%%% Siyu Xue -- Sep. 29, 2023
clear

%% Manually input the weighted avg velocities

A3 = [3.9 3.6 3.9 3.8 3.8 3.7];
A10 = [3.9 3.6 3.9 3.8 3.8 3.7];
A30 = [4.1 4 4.1 4 3.9 4];
A40 = [4.2 4.2 4.3 4.3 4.2 4];

B3 = [3.3 3.3 3.6];
B10 = [3.8 3.8 4];
B30 = [3.9 4.1 4];
B40 = [4.1 4.4 4.2];

M3 = [3.7 3.7 3.6];
M10 = [3.8 3.7 3.8];
M30 = [3.9 3.9 4];
M40 = [4.2 4.1 4.3];

O3 = [3.5 3.5];
O10 = [3.6 3.7];
O30 = [4 4.1];
O40 = [4.6 4.4];

U3 = [3.4 3.7 3.2];
U10 = [3.7 3.8 3.7];
U30 = [4.4 4 4.3];
U40 = [4.6 4.3 4.5];

%% Get the order
D3 = [A3 B3 M3 O3 U3];
D10 = [A10 B10 M10 O10 U10];
D30 = [A30 B30 M30 O30 U30];
D40 = [A40 B40 M40 O40 U40];

[D3_sort, D3_idx] = sort(D3, 'descend');
[D10_sort, D10_idx] = sort(D10, 'descend');
[D30_sort, D30_idx] = sort(D30, 'descend');
[D40_sort, D40_idx] = sort(D40, 'descend');

%% Plot all avg vels on the same figure
C1 = [0 0.447 0.741];
C2 = [0.851 0.325 0.098];
C3 = [0.929 0.694 0.125];
C4 = [0.467 0.674 0.188];

hold on

v1 = plot(D3_sort, 'Color', C1, LineWidth=3);
for i = 1:17

    idx = find(D3_idx == i);
  
    if i <= 6
        plot(idx, D3_sort(idx), "pentagram", 'MarkerFaceColor', C1, 'MarkerEdgeColor','k', 'MarkerSize',15)
        text(idx-0.2, D3_sort(idx)-0.05, num2str(i), 'Color', C1)
    elseif i <= 9
        plot(idx, D3_sort(idx), "o", 'MarkerFaceColor', C1, 'MarkerEdgeColor','k', 'MarkerSize',15)
        text(idx-0.2, D3_sort(idx)-0.05, num2str(i), 'Color', C1)
    elseif i <= 12
        plot(idx, D3_sort(idx), "square", 'MarkerFaceColor', C1, 'MarkerEdgeColor','k', 'MarkerSize',15)
        text(idx-0.2, D3_sort(idx)-0.05, num2str(i), 'Color', C1)
    elseif i <= 14
        plot(idx, D3_sort(idx), "^", 'MarkerFaceColor', C1, 'MarkerEdgeColor','k', 'MarkerSize',15)
        text(idx-0.2, D3_sort(idx)-0.05, num2str(i), 'Color', C1)
    else
        plot(idx, D3_sort(idx), "diamond", 'MarkerFaceColor', C1, 'MarkerEdgeColor','k', 'MarkerSize',15)
        text(idx-0.2, D3_sort(idx)-0.05, num2str(i), 'Color', C1)
    end
end

v2 = plot(D10_sort, 'Color', C2, LineWidth=3);
for i = 1:17

    idx = find(D10_idx == i);
  
    if i <= 6
        plot(idx, D10_sort(idx), "pentagram", 'MarkerFaceColor', C2, 'MarkerEdgeColor','k', 'MarkerSize',15)
        text(idx-0.05, D10_sort(idx)+0.05, num2str(i), 'Color', C2)
    elseif i <= 9
        plot(idx, D10_sort(idx), "o", 'MarkerFaceColor', C2, 'MarkerEdgeColor','k', 'MarkerSize',15)
        text(idx-0.05, D10_sort(idx)+0.05, num2str(i), 'Color', C2)
    elseif i <= 12
        plot(idx, D10_sort(idx), "square", 'MarkerFaceColor', C2, 'MarkerEdgeColor','k', 'MarkerSize',15)
        text(idx-0.05, D10_sort(idx)+0.05, num2str(i), 'Color', C2)
    elseif i <= 14
        plot(idx, D10_sort(idx), "^", 'MarkerFaceColor', C2, 'MarkerEdgeColor','k', 'MarkerSize',15)
        text(idx-0.05, D10_sort(idx)+0.05, num2str(i), 'Color', C2)
    else
        plot(idx, D10_sort(idx), "diamond", 'MarkerFaceColor', C2, 'MarkerEdgeColor','k', 'MarkerSize',15)
        text(idx-0.05, D10_sort(idx)+0.05, num2str(i), 'Color', C2)
    end
end

v3 = plot(D30_sort, 'Color', C3, LineWidth=3);
for i = 1:17

    idx = find(D30_idx == i);
  
    if i <= 6
        plot(idx, D30_sort(idx), "pentagram", 'MarkerFaceColor', C3, 'MarkerEdgeColor','k', 'MarkerSize',15)
        text(idx-0.2, D30_sort(idx)-0.06, num2str(i), 'Color', C3)
    elseif i <= 9
        plot(idx, D30_sort(idx), "o", 'MarkerFaceColor', C3, 'MarkerEdgeColor','k', 'MarkerSize',15)
        text(idx-0.2, D30_sort(idx)-0.06, num2str(i), 'Color', C3)
    elseif i <= 12
        plot(idx, D30_sort(idx), "square", 'MarkerFaceColor', C3, 'MarkerEdgeColor','k', 'MarkerSize',15)
        text(idx-0.2, D30_sort(idx)-0.06, num2str(i), 'Color', C3)
    elseif i <= 14
        plot(idx, D30_sort(idx), "^", 'MarkerFaceColor', C3, 'MarkerEdgeColor','k', 'MarkerSize',15)
        text(idx-0.2, D30_sort(idx)-0.06, num2str(i), 'Color', C3)
    else
        plot(idx, D30_sort(idx), "diamond", 'MarkerFaceColor', C3, 'MarkerEdgeColor','k', 'MarkerSize',15)
        text(idx-0.2, D30_sort(idx)-0.06, num2str(i), 'Color', C3)
    end
end

v4 = plot(D40_sort, 'Color', C4, LineWidth=3);
for i = 1:17

    idx = find(D3_idx == i);
  
    if i <= 6
        plot(idx, D40_sort(idx), "pentagram", 'MarkerFaceColor', C4, 'MarkerEdgeColor','k', 'MarkerSize',15)
        text(idx, D40_sort(idx)-0.06, num2str(i), 'Color', C4)
    elseif i <= 9
        plot(idx, D40_sort(idx), "o", 'MarkerFaceColor', C4, 'MarkerEdgeColor','k', 'MarkerSize',15)
        text(idx, D40_sort(idx)-0.06, num2str(i), 'Color', C4)
    elseif i <= 12
        plot(idx, D40_sort(idx), "square", 'MarkerFaceColor', C4, 'MarkerEdgeColor','k', 'MarkerSize',15)
        text(idx, D40_sort(idx)-0.06, num2str(i), 'Color', C4)
    elseif i <= 14
        plot(idx, D40_sort(idx), "^", 'MarkerFaceColor', C4, 'MarkerEdgeColor','k', 'MarkerSize',15)
        text(idx, D40_sort(idx)-0.06, num2str(i), 'Color', C4)
    else
        plot(idx, D40_sort(idx), "diamond", 'MarkerFaceColor', C4, 'MarkerEdgeColor','k', 'MarkerSize',15)
        text(idx, D40_sort(idx)-0.06, num2str(i), 'Color', C4)
    end
end

hold off
xlabel('Geological Region Index')
ylabel('Velocity (km/s)')
xticklabels([])
legend([v1, v2, v3, v4], 'Depth=3km','Depth=10km','Depth=30km','Depth=40km')
set(gcf,'color','w','position',[400,400,600,400]);
set(gca,'fontsize', 14);

% save the figure
figpath = '/Users/sxue3/Documents/BayMap_Figures/fig/figSX_CompareGeo.pdf';
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(gcf,figpath,'-painters', '-dpdf','-r0');

