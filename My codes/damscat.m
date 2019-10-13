

%% Scatter with damage gradient

function h = damscat(FileName, xcol, ycol, damageCol)

%% 
A = load(FileName);
A = sortrows(A, damageCol);
l = length(A);
dcol = A(:, damageCol);
CT=cbrewer('seq', 'YlOrBr', l);  % You have to download cbrewer file from Matlab File Exchange.


figure, 
scatter(A(:,xcol),A(:,ycol),47, dcol, 'filled');
colormap(CT);
% colormap(flipud(CT)); % To flip the colour map

h = gca;
a = colorbar('peer', h, 'Location', 'South','FontName','Times', ...
    'FontSize',15,'FontWeight','bold', 'TickLabels',{'Low','-','High'});

% Position and size of colourbar
b =  a.Position;
set(a,'Position',[(b(1)+b(3)/2) b(2)+b(4)/2 0.3 0.03]);


%% Giving title to colourbar

% cbtitle = get(a,'Title');
% title = '\gammaH2A.X';
% set(cbtitle ,'String',title,'FontName','Times','FontSize',17,'FontWeight','bold');

%% Nice Try though... makes graphs really really heavy.

% CT1 = cbrewer('seq', 'YlOrBr', floor(l/0.8));
% 
% for i = 1:l
% scatter(A(i,xcol),A(i,ycol),40,'MarkerEdgeColor',CT(i,:), 'MarkerFaceColor', CT1(i,:), 'LineWidth', 1.3);
% hold on;
% end