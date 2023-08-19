function out = ScatBar(cont_file, add_cont, damg_file, add_damg, col_num, varargin)


p = inputParser;

if nargin < 5
    error('MATLAB:narginchk:notEnoughInputs', 'Not enough input arguments.');
else
    addParameter(p, 'Scatter_on', 'n', ...
        @(t) (ischar(t) && ismember(t, {'n', 'y'})));
    addParameter(p, 'sig_test', 'n', ...
        @(t) (ischar(t) && ismember(t, {'ttest2','kstest2','n'})));    
end
p.KeepUnmatched = true;

parse(p, varargin{:});
Want = p.Results;
Want = Want.Scatter_on;
testor = p.Results;
testor = testor.sig_test;

% example: ScatBar(cont_file, add_cont, damg_file, add_damg, col_num, ...
% 'Scatter_on', 'y', 'sig_test', 'ttest2');

%% Getting control files


cd(add_cont)

if exist('G1.dat', 'file')...
   && exist('S.dat','file') && exist('G2_M.dat', 'file');
    
    fprintf('Good to go!\n\n\n');
    
    else
    
    fprintf('Lemme fix this.\n');
    CellCycleStage(cont_file);
    fprintf('Now you are good to go!\n\n\n');
    
end

u = load('G1.dat');     u_num = length(u);
v = load('S.dat');      v_num = length(v);
w = load('G2_M.dat');   w_num = length(w);


%% Getting other files


cd(add_damg)

if exist('G1.dat', 'file')...
   && exist('S.dat','file') && exist('G2_M.dat', 'file');
    
    fprintf('Good to go!\n\n\n');
    
    else
    
    fprintf('Lemme fix this.\n');
    CellCycleStage(damg_file);
    fprintf('Now you are good to go!\n\n\n');
    
end

x = load('G1.dat');     x_num = length(x);
y = load('S.dat');      y_num = length(y);
z = load('G2_M.dat');   z_num = length(z);

%% Plotting bar graph
 
 a = mean(u); b = mean(v); c = mean(w);
 d = std(u)/sqrt(length(u)); e = std(v)/sqrt(length(v)); f = std(w)/sqrt(length(w));
 g = mean(x); h = mean(y); i = mean(z);
 j = std(x)/sqrt(length(x)); k = std(y)/sqrt(length(y)); l = std(z)/sqrt(length(z));

 %%
 p = [a(:,col_num) g(:,col_num); b(:,col_num) h(:,col_num); c(:,col_num) i(:,col_num)];
 q = [d(:,col_num) j(:,col_num); e(:,col_num) k(:,col_num); f(:,col_num) l(:,col_num)];
axes1 = axes('Parent',figure);
 xOrder = 1:size(p,1);
 [~, nCols] = size(p);
% CT=cbrewer('seq', 'YlOrBr', 3);  % can change it to RdYlGn/RdYlBu try it.
CT = [250/255,240/255,0/255; 200/255,115/255,0/255];
bwe = bar(p, 'LineWidth', 3);
set(gca,'XTick',[1 2 3],'XTickLabel',{'G1','S','G2/M'});
colormap(CT);
 hold on;
 hErrorbar = zeros(1,nCols);
 
 set(axes1,'FontName','Times','FontSize',37,'FontWeight','bold')

 %% For scatter on bars
 
Scatter_on = Want;
    switch Scatter_on
        case 'y'
        for i = 1:nCols
      
            bl = bwe(i).BarWidth/4;
            xdata =  bwe(i).XData + [bwe(i).XOffset];
      
      
            if i == 1
          
                g1_x = (xdata(1,1)-bl/2):(bl/(u_num-1)):(xdata(1,1)+bl/2);
                s_x = (xdata(1,2)-bl/2):(bl/(v_num-1)):(xdata(1,2)+bl/2);
                g2_x = (xdata(1,3)-bl/2):(bl/(w_num-1)):(xdata(1,3)+bl/2);
                s1 = scatter('v6', g1_x,u(:,col_num),5,'MarkerEdgeColor',[.2*i .5 .7], 'MarkerFaceColor',[1 .0 .0], 'LineWidth',.5);
                s2 = scatter('v6', s_x,v(:,col_num),5,'MarkerEdgeColor',[.3*i .6 .8], 'MarkerFaceColor',[1 .0 .0], 'LineWidth',.5);
                s3 = scatter('v6', g2_x,w(:,col_num),5,'MarkerEdgeColor',[.4*i .7 .9], 'MarkerFaceColor',[1 .0 .0], 'LineWidth',.5);
                s1_f = get(s1, 'children');
                s2_f = get(s2, 'children');
                s3_f = get(s3, 'children');
                set(s1_f, 'FaceAlpha', 0.3, 'EdgeAlpha', 0.2);
                set(s2_f, 'FaceAlpha', 0.3, 'EdgeAlpha', 0.2);
                set(s3_f, 'FaceAlpha', 0.3, 'EdgeAlpha', 0.2);
            else
               
                g1_x = (xdata(1,1)-bl/2):(bl/(x_num-1)):(xdata(1,1)+bl/2);
                s_x = (xdata(1,2)-bl/2):(bl/(y_num-1)):(xdata(1,2)+bl/2);
                g2_x = (xdata(1,3)-bl/2):(bl/(z_num-1)):(xdata(1,3)+bl/2);
                s1 = scatter('v6', g1_x,x(:,col_num),5,'MarkerEdgeColor',[.2*i .5 .7], 'MarkerFaceColor',[.2*i .4 .6], 'LineWidth',.5);
                s2 = scatter('v6', s_x,y(:,col_num),5,'MarkerEdgeColor',[.3*i .6 .8], 'MarkerFaceColor',[.3*i .5 .7], 'LineWidth',.5);
                s3 = scatter('v6', g2_x,z(:,col_num),5,'MarkerEdgeColor',[.4*i .7 .9], 'MarkerFaceColor',[.4*i .6 .8], 'LineWidth',.5);
                s1_f = get(s1, 'children');
                s2_f = get(s2, 'children');
                s3_f = get(s3, 'children');
                set(s1_f, 'FaceAlpha', 0.2, 'EdgeAlpha', 0.2);
                set(s2_f, 'FaceAlpha', 0.2, 'EdgeAlpha', 0.2);
                set(s3_f, 'FaceAlpha', 0.2, 'EdgeAlpha', 0.2);
          
            end
      
            hErrorbar(i) = errorbar(mean(xdata,1), p(xOrder,i), q(xOrder,i), q(xOrder, i), '.k');
            set(hErrorbar(i), 'marker', 'none', 'LineWidth', 3)
        end
        
        otherwise
            
            for i = 1:nCols
                xdata =  bwe(i).XData + [bwe(i).XOffset];
                hErrorbar(i) = errorbar(mean(xdata,1), p(xOrder,i), q(xOrder,i), q(xOrder, i), '.k');
                set(hErrorbar(i), 'marker', 'none', 'LineWidth', 3)
            end
    end
 
 %% Significance test 
 
 sig_test = testor;
 sigwant = isequal(testor, 'n');
 sigwant = double(sigwant);
 
 if sigwant == 0
     
 pv = zeros(1,3); 
 switch sig_test
         
     case 'ttest2'
        [~, pv(1)] = ttest2(u(:,col_num), x(:,col_num));
        [~, pv(2)] = ttest2(v(:,col_num), y(:,col_num));
        [~, pv(3)] = ttest2(w(:,col_num), z(:,col_num));
         
     case 'kstest2'
         [~, pv(1)] = kstest2(u(:,col_num), x(:,col_num));
         [~, pv(2)] = kstest2(v(:,col_num), y(:,col_num));
         [~, pv(3)] = kstest2(w(:,col_num), z(:,col_num));
     otherwise
 end
     
 for i = 1:3
     
     if p(i,1) > p(i,2)
         bn = 1;
     else
         bn = 2;
     end
      A = max(p(i,1), p(i,2));
      B = min(p(i,1), p(i,2));
      C = max([q(1,1), q(1,2), q(2,1), q(2,2), q(3,1), q(3,2)]);
      C = min(C, (A-B));
     
     if pv(i) >= 0.05
         
     elseif pv(i) >= 0.00005
            line([(bwe(1).XData(i) - bwe(bn).XOffset) (bwe(1).XData(i) - bwe(bn).XOffset)], [(B+2*C) (A+3*C)],'LineWidth',3, 'color', 'k');
            line([(bwe(1).XData(i) - bwe(bn).XOffset) (bwe(1).XData(i) + bwe(bn).XOffset)], [(A+3*C) (A+3*C)],'LineWidth',3, 'color', 'k');
            line([(bwe(1).XData(i) + bwe(bn).XOffset) (bwe(1).XData(i) + bwe(bn).XOffset)], [(A+3*C) (A+2*C)],'LineWidth',3, 'color', 'k');
            text((bwe(1).XData(i) - 1.5*bwe(2).XOffset), (A+4*C), '\bf *', 'FontSize',19, 'color', 'k');
     elseif pv(i) >= 0.0000005
            line([(bwe(1).XData(i) - bwe(bn).XOffset) (bwe(1).XData(i) - bwe(bn).XOffset)], [(B+2*C) (A+3*C)],'LineWidth',3, 'color', 'k');
            line([(bwe(1).XData(i) - bwe(bn).XOffset) (bwe(1).XData(i) + bwe(bn).XOffset)], [(A+3*C) (A+3*C)],'LineWidth',3, 'color', 'k');
            line([(bwe(1).XData(i) + bwe(bn).XOffset) (bwe(1).XData(i) + bwe(bn).XOffset)], [(A+3*C) (A+2*C)],'LineWidth',3, 'color', 'k');
            text((bwe(1).XData(i) - 1.5*bwe(2).XOffset), (A+4*C), '\bf **', 'FontSize',19, 'color', 'k');
     else
            line([(bwe(1).XData(i) - bwe(bn).XOffset) (bwe(1).XData(i) - bwe(bn).XOffset)], [(B+2*C) (A+3*C)],'LineWidth',3, 'color', 'k');
            line([(bwe(1).XData(i) - bwe(bn).XOffset) (bwe(1).XData(i) + bwe(bn).XOffset)], [(A+3*C) (A+3*C)],'LineWidth',3, 'color', 'k');
            line([(bwe(1).XData(i) + bwe(bn).XOffset) (bwe(1).XData(i) + bwe(bn).XOffset)], [(A+3*C) (A+2*C)],'LineWidth',3, 'color', 'k');
            text((bwe(1).XData(i)  - 1.5*bwe(2).XOffset), (A+4*C), '\bf ***', 'FontSize',19, 'color', 'k');
                
       
     end
 end
 
 else
     
 end
hold off


end