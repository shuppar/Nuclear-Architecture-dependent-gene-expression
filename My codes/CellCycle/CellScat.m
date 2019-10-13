%% Fitting bimodality to the DNA Content.

% Number of bins is calculated by Scott's rule.

function h = CellScat(filename, xcol, ycol, varargin)


%%


p = inputParser;

if nargin < 3
    error('MATLAB:narginchk:notEnoughInputs', 'Not enough input arguments.');
else
    addParameter(p, 'Normalize', 'y', ...
        @(t) (ischar(t) && ismember(t, {'n', 'y'})));  
end
p.KeepUnmatched = true;

parse(p, varargin{:});
Norm = p.Results;
Norm = Norm.Normalize;



%%

if exist('G1_peak.dat', 'file') && exist('G1.dat', 'file')...
        && exist('S.dat','file') && exist('G2_M.dat', 'file');
%         && exist('BeyondG2_M.dat', 'file') && exist('subG1.dat', 'file');
    
    fprintf('Good to go!\n\n\n');
    
    
else
    
    fprintf('Lemme fix this.\n');
    CellCycleStage(filename);
    fprintf('Now you are good to go!\n\n\n');
    
end

    A4 = load(filename);
    G = load('G1_peak.dat'); G = G(:,1);
    C = load('S.dat'); B = load('G1.dat'); D = load('G2_M.dat');
    B(:,1) = B(:,1)/G; C(:,1) = C(:,1)/G; D(:,1) = D(:,1)/G; A4(:,1) = A4(:,1)/G;    % Normalization wrt G1 peak added on September 11, 2017
    
    
    if exist('BeyondG2_M.dat', 'file') && exist('subG1.dat', 'file')
        A = load('subG1.dat'); E = load('BeyondG2_M.dat');      % Normalization wrt G1 peak added on September 11, 2017
        A(:,1) = A(:,1)/G; E(:,1) = E(:,1)/G;        
    else
        
    end
   
%% Normalization yes/no?
    
switch Norm
    
    case 'y'

        if xcol == 1 && ycol~= 1
            meanx = 1; meany = mean(A4(:,ycol));
        elseif ycol == 1 && xcol ~=1
            meany = 1; meanx = mean(A4(:,xcol));
        elseif xcol == 1 && ycol == 1
            meanx = 1; meany = 1;
        elseif xcol ~= 1 && ycol ~= 1
            meanx = mean(A4(:,xcol)); meany = mean(A4(:,ycol));
        end
    
    
        B(:,xcol) = B(:,xcol)/meanx; B(:,ycol) = B(:,ycol)/meany;
        C(:,xcol) = C(:,xcol)/meanx; C(:,ycol) = C(:,ycol)/meany;
        D(:,xcol) = D(:,xcol)/meanx; D(:,ycol) = D(:,ycol)/meany;
    
            if exist('BeyondG2_M.dat', 'file') && exist('subG1.dat', 'file')
                A(:,xcol) = A(:,xcol)/meanx; A(:,ycol) = A(:,ycol)/meany;
                E(:,xcol) = E(:,xcol)/meanx; E(:,ycol) = E(:,ycol)/meany;
            else
        
            end
            
    otherwise
end
    
    
    %% Plotting
    
    

    axes1 = axes('Parent',figure);
    
    scatter(B(:,xcol),B(:,ycol),40,'MarkerEdgeColor',[254,216,0]/255, 'MarkerFaceColor',[254,246,0]/255, 'LineWidth',1.5);
    h = gca;
    hold on;
    scatter(C(:,xcol),C(:,ycol),40,'MarkerEdgeColor',[235,110,0]/255, 'MarkerFaceColor', [250,155,0]/255, 'LineWidth',1.5);
    scatter(D(:,xcol),D(:,ycol),40,'MarkerEdgeColor',[80,17,6]/255, 'MarkerFaceColor',[150,40,20]/255, 'LineWidth',1.5);
    legend('G1', 'S', 'G2/M');
    if exist('BeyondG2_M.dat', 'file') && exist('subG1.dat', 'file')
        scatter(A(:,xcol),A(:,ycol),40,'MarkerEdgeColor',[254,216,0]/255, 'MarkerFaceColor',[254,246,0]/255, 'LineWidth',1.5);
        scatter(E(:,xcol),E(:,ycol),40,'MarkerEdgeColor',[80,17,6]/255, 'MarkerFaceColor',[150,40,20]/255, 'LineWidth',1.5);
    else
    end
    
    hold off;
    
xlabel({'(a.u.)'}); ylabel({'(a.u.)'});
% legend(axes1,'show', 'Orientation', 'Horizontal');
box(axes1,'on');
set(axes1,'FontName','Times','FontSize',37,'FontWeight','bold');
title('Title', 'FontSize', 37, 'FontWeight', 'bold', 'FontName', 'Times');
    
end
