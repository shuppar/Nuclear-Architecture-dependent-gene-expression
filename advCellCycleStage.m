%% Shuppar script for segregating cells in different stages of cell cycle

% Requires EdU staining

function [] = advCellCycleStage(FileName, EdU_column_no, varargin)

%%

p = inputParser;

if nargin < 2
    error('MATLAB:narginchk:notEnoughInputs', 'Not enough input arguments.');
else
    addParameter(p, 'thresh', 1.3', ...
        @(t) (isnumeric(t) && isreal(t)));  
end
p.KeepUnmatched = true;

parse(p, varargin{:});thresh = p.Results;
t = thresh.thresh;


%% 

num = EdU_column_no;

    if exist('S.dat', 'file')
    
            A = load('S.dat');
    else
            fprintf('Lemme fix this.\n');
            CellCycleStage(FileName);
            A = load('S.dat');
    end
    
%     AS = sort(A(:,num), 'descend'); las = floor(length(AS)/2);
%     AS = mean(AS(1:las)); % mean of first two quartiles.
    AS = mean(A(:,num));
    SS = std(A(:, num)); SD = t*SS;
    peaks = load('G1_peak.dat');
    g1p = peaks(1,1);
    s1p = peaks(1,3);
    g2p = peaks(1,2);
    s2p = peaks(1,4);
    
G1 = load('G1.dat'); G2 = load('G2_M.dat');

for i = 1:size(G1,1)
   
    if (G1(i,num) < (AS-SD)) % Can change the multiplication factor (F) in f*SS
        
       f = fopen('newG1.dat','a');  
       fprintf(f,'%d\t', G1(i,:));
       fprintf(f,'\n');

       
    else
        f = fopen('newS.dat','a');  
        fprintf(f,'%d\t', G1(i,:));
	    fprintf(f,'\n');
        
               
       f = fopen('ES.dat','a');  
       fprintf(f,'%d\t', G1(i,:));
       fprintf(f,'\n');
        
    end
    
    
end


for i = 1:size(G2,1)
   
    if (G2(i,num) < (AS-SD))
        
       f = fopen('newG2_M.dat','a');  
       fprintf(f,'%d\t', G2(i,:));
       fprintf(f,'\n');

    else
        f = fopen('newS.dat','a');  
        fprintf(f,'%d\t', G2(i,:));
	    fprintf(f,'\n');
        
               
       f = fopen('LS.dat','a');  
       fprintf(f,'%d\t', G2(i,:));
       fprintf(f,'\n');
        
    end
    
    
end


for i = 1:size(A,1)
   
    if (A(i,num) < (AS-SD))
        
        if (abs(A(i,1)-g1p) <= abs(A(i,1)-g2p))
        
        f = fopen('newG1.dat','a');  
        fprintf(f,'%d\t', A(i,:));
	    fprintf(f,'\n');
            
             
        else
            
        f = fopen('newG2_M.dat','a');  
        fprintf(f,'%d\t', A(i,:));
	    fprintf(f,'\n');
            
        end
       
    else
        f = fopen('newS.dat','a');  
        fprintf(f,'%d\t', A(i,:));
	    fprintf(f,'\n');
        
    end
    
    
end

if exist('subG1.dat', 'file')
    
    lG1 = load('subG1.dat');
        for i = 1:size(lG1,1)
   
            if (lG1(i,num) < (AS-SD)) % Can change the multiplication factor (F) in f*SS
        
                    f = fopen('newsubG1.dat','a');  
                    fprintf(f,'%d\t', lG1(i,:));
                    fprintf(f,'\n');
                           
                else
                    f = fopen('newS.dat','a');  
                    fprintf(f,'%d\t', lG1(i,:));
                    fprintf(f,'\n');
                    
                    f = fopen('ES.dat','a');  
                    fprintf(f,'%d\t', lG1(i,:));
                    fprintf(f,'\n');
        
            end
    
    
        end
        
        if exist('newsubG1.dat', 'file')
        movefile('newsubG1.dat', 'subG1.dat');
        end
        
        
else

end

if exist('BeyondG2_M.dat', 'file')
    
    rG2 = load('BeyondG2_M.dat');
        for i = 1:size(rG2,1)
   
            if (rG2(i,num) < (AS-SD)) % Can change the multiplication factor (F) in f*SS
        
                    f = fopen('newBG2_M.dat','a');  
                    fprintf(f,'%d\t', rG2(i,:));
                    fprintf(f,'\n');
                           
                else
                    f = fopen('newS.dat','a');  
                    fprintf(f,'%d\t', rG2(i,:));
                    fprintf(f,'\n');
                    
                    f = fopen('LS.dat','a');  
                    fprintf(f,'%d\t', rG2(i,:));
                    fprintf(f,'\n');
        
            end
    
    
        end
        
        if exist('newBG2_M.dat', 'file')
        movefile('newBG2_M.dat', 'BeyondG2_M.dat');
        end
        
else

end

NS = load('newS.dat');

for i = 1:size(NS,1)
   
    if (NS(i,1) > (g1p+2*s1p) && NS(i,1) < (g2p-0.8*s1p))
        
       f = fopen('MS.dat','a');  
       fprintf(f,'%d\t', NS(i,:));
       fprintf(f,'\n');
       
    else
        
    end
    
    
end


movefile('newG1.dat', 'G1.dat');
movefile('newS.dat', 'S.dat');
movefile('newG2_M.dat', 'G2_M.dat');

end
