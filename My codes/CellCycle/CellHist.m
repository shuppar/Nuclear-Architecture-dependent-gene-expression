%% Fitting bimodality to the DNA Content.

% Number of bins is calculated by Scott's rule.

function h = CellHist(filename)


%%
if exist('G1_peak.dat', 'file') && exist('G1.dat', 'file')...
        && exist('S.dat','file') && exist('G2_M.dat', 'file');
%         && exist('BeyondG2_M.dat', 'file') && exist('subG1.dat', 'file');
    
    fprintf('Good to go!\n');
    
    
else
   
    CellCycleStage(filename);
    fprintf('Now you are good to go!\n');
    
end


    
    C = load('S.dat');
    B = load('G1.dat'); D = load('G2_M.dat'); 
    F = load(filename); G = load('G1_peak.dat');
    B(:,1) = B(:,1)/G(:,1); C(:,1) = C(:,1)/G(:,1);     % Normalization wrt G1 peak added on September 11, 2017
    D(:,1) = D(:,1)/G(:,1); F(:,1) = F(:,1)/G(:,1);     % Normalization wrt G1 peak added on September 11, 2017
    
%%
    
    bl = (3.5*std(F(:,1)))/((length(F))^(1/3));
 
    h2 = (max(B(:,1)) - min(B(:,1)))/bl; h2 = round(h2); 
    h3 = (max(C(:,1)) - min(C(:,1)))/bl; h3 = round(h3); 
    h4 = (max(D(:,1)) - min(D(:,1)))/bl; h4 = round(h4);  
    h6 = (max(F(:,1)) - min(F(:,1)))/bl; h6 = round(h6); 
    
    if exist('BeyondG2_M.dat', 'file') && exist('subG1.dat', 'file')
        A = load('subG1.dat'); E = load('BeyondG2_M.dat');
        A(:,1) = A(:,1)/G(:,1); E(:,1) = E(:,1)/G(:,1);          % Normalization wrt G1 peak added on September 11, 2017
        h1 = (max(A(:,1)) - min(A(:,1)))/bl; h1 = round(h1);
        h5 = (max(E(:,1)) - min(E(:,1)))/bl; h5 = round(h5);
    else
        
    end
    
    G(:,3) = G(:,3)/G(:,1); G(:,4) = G(:,4)/G(:,1); % Normalization wrt G1 peak added on September 11, 2017
    G(:,2) = G(:,2)/G(:,1); G(:,1) = G(:,1)/G(:,1); % Normalization wrt G1 peak added on September 11, 2017
       
    u = F(:,1);
    v1 = (G(:,1) - bl/2); w1 = (G(:,1) + bl/2);
    s1 = sum(u>v1 & u<w1);
    t1 = s1/G(:,5);
    v2 = (G(:,2) - bl/2); w2 = (G(:,2) + bl/2);
    s2 = sum(u>v2 & u<w2); 
    t2 = s2/(1-G(:,5));
    ugrid = linspace(0.8*min(u),2*max(u),200);
    bimode = @(x,p,mu1,mu2,sigma1,sigma2) (sqrt(2*pi)*sigma1)*t1*p*normpdf(x,mu1,sigma1) + (sqrt(2*pi)*sigma2)*t2*(1-p)*normpdf(x,mu2,sigma2);
    pdfgrid = bimode(ugrid,G(:,5),G(:,1),G(:,2),G(:,3), G(:,4));
    
    %% Separate fitting... May 22, 2017
    
    if exist('G2_peak.dat', 'file')
    G2 = load('G2_peak.dat');
    s2gauss = @(x,mu1,sigma1,s) (s*(sqrt(2*pi)*sigma1)*normpdf(x,mu1,sigma1));
    sugrid1 = linspace(0.8*min(D(:,1)), 2*max(D(:,1)), 200);
    spdfgrid1 = s2gauss(sugrid1, G2(:,1), G2(:,2), s2);
    sugrid2 = linspace(0.8*min(B(:,1)), 2*max(B(:,1)), 200);
    spdfgrid2 = s2gauss(sugrid2, G(:,1), G(:,3), s1);
    
    else
        
        SCellFit(filename);
        G2 = load('G2_peak.dat');
        s2gauss = @(x,mu1,sigma1,s) (s*(sqrt(2*pi)*sigma1)*normpdf(x,mu1,sigma1));
        sugrid1 = linspace(0.8*min(D(:,1)), 2*max(D(:,1)), 200);
        spdfgrid1 = s2gauss(sugrid1, G2(:,1), G2(:,2), s2);
        sugrid2 = linspace(0.8*min(B(:,1)), 2*max(B(:,1)), 200);
        spdfgrid2 = s2gauss(sugrid2, G(:,1), G(:,3), s1);
        
    end
    
    %% Cell Cycle Plot
    
    x1 = mean(B(:,1)); y1 = (length(B(:,1))/4);
    x2 = mean(D(:,1)); y2 = (length(D(:,1))/4);
%     cv1 = std(B(:,1))/x1; cv2 = std(D(:,1))/x2;  % From arbitrary definition of G2- phase cell.
    

    cv1 = G(:,3)/G(:,1); cv2 = G(:,4)/(G(:,2));       % From bi-modal fit.
    
    figure, histogram(B(:,1), max(h2,1)); xlabel('DNA Content (a.u.)', 'FontSize', 37, 'FontName', 'Times', 'FontWeight', 'Bold'); ylabel('Count', 'FontSize', 37, 'FontName', 'Times', 'FontWeight', 'Bold');
    hold on
    histogram(C(:,1), max(h3,1));
    histogram(D(:,1), max(h4,1));
    
    if exist('A', 'var') && exist('E', 'var')
        histogram(A(:,1), max(h1,1));
        histogram(E(:,1), max(h5,1));
    else
    end
    plot(ugrid,pdfgrid,'-');  

    text(x1, y1, sprintf('%.2f', cv1));
    text(x2, y2, sprintf('%.2f', cv2));    
    
        env1 = plot(ugrid,pdfgrid,'-');
        if exist('G2_peak.dat', 'file')
        
        env2 = plot(sugrid1, spdfgrid1, '-');
        plot(sugrid2, spdfgrid2, '-');
        legend([env1 env2],{'Bimodal','Separate'}); 
        else
        end
    
    hold off;    

    %% Histogram Plot
    
    figure, histogram(u, 'LineWidth', 3, 'EdgeColor', [235,110,0]/255, 'FaceColor', [250,155,0]/255, 'FaceAlpha', 0.6, 'EdgeAlpha', 0.6);
    xlabel('DNA Content (a.u.)', 'FontSize', 37, 'FontName', 'Times', 'FontWeight', 'Bold'); ylabel('Count', 'FontSize', 37,'FontName', 'Times', 'FontWeight', 'Bold');
    hold on;
    env1 = plot(ugrid,pdfgrid,'-');
    text(x1, y1, sprintf('%.2f', cv1));
    text(x2, y2, sprintf('%.2f', cv2)); 
    
    if exist('G2_peak.dat', 'file')
        
        env2 = plot(sugrid1, spdfgrid1, '-');
        plot(sugrid2, spdfgrid2, '-');
        legend([env1 env2],{'Bimodal','Separate'}); 
    else
    end
    
    hold off;
    
end
