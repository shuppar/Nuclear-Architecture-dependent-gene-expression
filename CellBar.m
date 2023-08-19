%% Fitting bimodality to the DNA Content.

% Number of bins is calculated by Scott's rule.

function h = CellBar(filename, xcol, varargin)


%%


p = inputParser;

if nargin < 2
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

        if xcol == 1
            meanx = 1;
        else
            meanx = mean(A4(:,xcol));
        end
    
    
        B(:,xcol) = B(:,xcol)/meanx; 
        C(:,xcol) = C(:,xcol)/meanx; 
        D(:,xcol) = D(:,xcol)/meanx; 
    
%             if exist('BeyondG2_M.dat', 'file') && exist('subG1.dat', 'file')
%                 A(:,xcol) = A(:,xcol)/meanx; A(:,ycol) = A(:,ycol)/meany;
%                 E(:,xcol) = E(:,xcol)/meanx; E(:,ycol) = E(:,ycol)/meany;
%             else
%         
%             end
            
    otherwise
end
    
    
    %% Plotting
    
    

    axes1 = axes('Parent',figure);
    
    histogram(B(:,xcol), 'FaceAlpha', 0.95, 'EdgeColor',[254,216,0]/255, 'FaceColor', [254,246,0]/255, 'LineWidth', 3);
    h = gca;
    hold on;
    histogram(C(:,xcol), 'FaceAlpha', 0.95, 'EdgeColor',[235,110,0]/255, 'FaceColor', [250,155,0]/255, 'LineWidth', 3);
    histogram(D(:,xcol), 'FaceAlpha', 0.95, 'EdgeColor',[80,17,6]/255, 'FaceColor', [150,40,20]/255, 'LineWidth', 3);
%     legend('G1', 'S', 'G2/M');    
    hold off;
    
xlabel({'(a.u.)'}); ylabel({'# of Cells'});
% legend(axes1,'show', 'Orientation', 'Horizontal');
box(axes1,'on');
set(axes1,'FontName','Times','FontSize',37,'FontWeight','bold');
title('Title', 'FontSize', 37, 'FontWeight', 'bold', 'FontName', 'Times');
    
end
