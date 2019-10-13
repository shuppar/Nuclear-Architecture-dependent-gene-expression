% Overlap statistics

% Finds percentage of overlap between two mainly punctated/spotty images.
% Works with any two binary images. But better if you give cropped image of a region
% you want to find overlap in.

function [num1, num1o, num1no, POBI1, num2, num2o, num2no, POBI2] = OvLap(BI1, BI2, MinOP)

% MinOA = minimum overlap area smaller than which will not be counted for
% the statistics.
% BI1 = binary Image 1, similarly for BI2.
% POBI1 = percentage overlap for BI1, similarly for POBI2.

I1 = BI1; I2 = BI2; ma = MinOP;

%% Some funny thing

% bgmark1 = -bwdist(I1); bgmark2 = -bwdist(I2);
% bgmark1(~I1) = Inf; bgmark2(~I2) = Inf;
% bgmark1 = watershed(bgmark1); bgmark2 = watershed(bgmark2);
% bgmark1(~I1) = 0; bgmark2(~I2) = 0;
% I1 = bgmark1; I2 = bgmark2; clear bgmark1; clear bgmark2

%%

[nrows, ncols] = size(I1);
I1 = logical(I1); I2 = logical(I2);
[~, num1] = bwlabeln(I1);
[~, num2] = bwlabeln(I2);
I1p = double(I1);
I2p = double(I2);

Io = I1p + I2p;
Io(Io~=2) = 0; Io = logical(Io);

R1 = regionprops(I1, 'Area', 'BoundingBox');
R2 = regionprops(I2, 'Area', 'BoundingBox');

POBI1 = 0; num1o = 0; num1no = 0;
POBI2 = 0; num2o = 0; num2no = 0;
%%

if num1 == 0

else

    for i = 1:num1
        xc = ceil(R1(i).BoundingBox);
        yc = shuppar(xc(2)-5); xw = dhuppar(xc(3)+5, ncols);
        yw = dhuppar(xc(4)+5, nrows); xc = (xc(1)-5);
        nL1 = Io(yc:yc+yw, xc:xc+xw, :); % figure, imshow(nL1);
        tr1 = regionprops(nL1, 'Area');
        [~, nt] = bwlabeln(nL1);
        ta = 0;
        for j = 1:nt
            ta = ta + tr1(j).Area;
        end

        if ta < R1(i).Area*ma
        else
            num1o = num1o+1;
        end
        
    end
    
    num1no = num1 - num1o;
    POBI1 = (num1o*100)/num1;
end


%%
if num2 == 0

else

    for i = 1:num2
        xc = ceil(R2(i).BoundingBox);
        yc = shuppar(xc(2)-5); xw = dhuppar(xc(3)+5, ncols);
        yw = dhuppar(xc(4)+5, nrows); xc = (xc(1)-5);
        nL1 = Io(yc:yc+yw, xc:xc+xw, :);
        tr1 = regionprops(nL1, 'Area'); % figure, imshow(nL1);
        [~, nt] = bwlabeln(nL1);
        ta = 0;
        for j = 1:nt
            ta = ta + tr1(j).Area;
        end
        
        if ta < R2(i).Area*ma
        else
            num2o = num2o+1;
        end
        
    end
    
    num2no = num2 - num2o;
    POBI2 = (num2o*100)/num2;
end

