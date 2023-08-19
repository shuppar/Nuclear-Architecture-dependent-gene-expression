%% A few things to note (Important)...

% Do flat field correction...
% it's simple read more on wikipedia or 
% @ nic.ucsf.edu/resources/how-to-acquire-flat-field-correction-images/
% So you might need to chage the disk size in strel in the second section.
% Calculate the minimum and maximum nuclear sizes and make changes
% accordingly in the section devoted to it.


%% Shuppar Script for Nuclear Intensity Plots.

close all;
clear all;
clc;

fmain=sprintf('18');
fileFolder = fullfile('18');
dirDAPI = dir(fullfile(fileFolder, '*C0001.tif')); 
dirEU = dir(fullfile(fileFolder, '*C0002.tif'));
dirHist = dir(fullfile(fileFolder, '*C0003.tif'));
dirRPA = dir(fullfile(fileFolder, '*C0004.tif'));
[DAPI_S, DAPI] = AvgNStack(fileFolder, dirDAPI);
[EU_S, ~] = AvgNStack(fileFolder, dirEU);
[hist_S, ~] = AvgNStack(fileFolder, dirHist);
[RPA_S, ~] = AvgNStack(fileFolder, dirRPA);

[xs, ys, zs] = size(DAPI_S); tot = xs*ys; clear xs; clear ys;
[nrows, ncols] = size(DAPI);

%% Nonuniform illumination correction. Function saved in home/MatlabCode. Copy of Image Analyst's code. (August 2, 2016)

BlankImage = imread('BlankImage.tif');

DAPI_S = BackgroundCorrect(DAPI_S, BlankImage);
EU_S = BackgroundCorrect(EU_S, BlankImage);
hist_S = BackgroundCorrect(hist_S, BlankImage);
RPA_S = BackgroundCorrect(RPA_S, BlankImage);

DAPI = uint16(mean(DAPI_S, 3));
B = uint16(mean(EU_S, 3));
C = uint16(mean(RPA_S, 3));
D = uint16(mean(hist_S, 3));
Blank = zeros(size(DAPI));
% phal = phal + C/1.5;

%% Rolling ball subtraction

% DAPI = imtophat(DAPI, strel('disk', 121));
% B = imtophat(B, strel('disk', 121));
% C = imtophat(C, strel('disk', 121));
% A = DAPI;

%% Subtracting Background

if exist('bg.dat', 'file')
    
    bagr = load('bg.dat');
    DB = bagr(1); HB = bagr(2); EB = bagr(3); RB = bagr(4);
    
else

M = imdilate(im2bw(DAPI, graythresh(DAPI)), strel('disk', 51));
M = ~M;
[~, bgo] = bwlabel(M);
% bgo
bgrnd = regionprops(M, 'Area');
Area = 0;

for i = 1:bgo
    Area = Area + bgrnd(i).Area;
end

M = uint16(M);
DB = max(DAPI_S, [], 3).*(M);
EB = max(EU_S, [], 3).*(M);
RB = max(RPA_S, [], 3).*(M);
HB = max(hist_S, [], 3).*(M);

DB = (mean(DB(:))*tot)/(Area); HB = (mean(HB(:))*tot)/(Area);
EB = (mean(EB(:))*tot)/(Area); RB = (mean(RB(:))*tot)/(Area);
% CB = (mean(CB(:))*tot)/(Area);

f = fopen('bg.dat', 'w');  
fprintf(f,'%f\t%f\t%f\t%f\n', DB, HB, EB, RB);


end

DAPI = DAPI - DB; D = D - HB; C = C -RB; B = B - EB; 
DAPI(DAPI<1) = 0; B(B<1) = 0; D(D<1) = 0; C(C<1) = 0;
% C = C - CB;
A = DAPI;

%% Nuclear Mask

[L, ~] = NMask(A, 11, 8500, 25000, 0.25, 1.4, 0.25, 1.4);  % NMask(Image, DiskRadius, minArea, maxArea, minCircularity, maxCircularity, minRoundness, maxRoundness)
L = mimclearborder(L, 10); % Just for those cases where cells are to be segmented.
[~, num] = bwlabel(L);
% Circularity is bigger than Roundness (at least for an ellipse) but ...
% the way they are defined are always smaller than 1.
% Ellipse; because most nuclei are elliptic.

% num = numel(cm1);
fprintf('\n\nThe number of cells detected is: %d\n\n', num);
NumLab = regionprops(L, 'Centroid', 'PixelIdxList', 'BoundingBox', 'Area');
% visMask = visMask(L);


%% Getting mean intensity of DAPI, H2AX and Ki67 for the detected cells.

dna = regionprops(L, DAPI, 'Area', 'MeanIntensity');
EU = regionprops(L, B, 'MeanIntensity');
RPA = regionprops(L, C, 'MeanIntensity');
% rp2 = regionprops(L, C, 'Area', 'MeanIntensity');

hist = regionprops(L, D, 'Area','MeanIntensity');


%% Colouring Cells in different stages of cell cycle in different colours

cd Data
peak1 = load('G1_peak.dat');
g1 = peak1(1,1);
s1 = peak1(1,3);
g2 = peak1(1,2);
s2 = peak1(1,4);
A = load('S.dat');
AS = mean(A(:,6));
SS = std(A(:, 6)); SD = 1.3*SS;
Red = zeros(size(L));
Blue = Red; Green = Red;

for i = 1:num

if (EU(i).MeanIntensity*dna(i).Area) < (AS-SD)

        elseif (dna(i).MeanIntensity*NumLab(i).Area < (g1+2.0*s1))
                 Red(NumLab(i).PixelIdxList) = 1;
                 Green(NumLab(i).PixelIdxList) = 0;
                 Blue(NumLab(i).PixelIdxList) = 0;
                 
                 
             elseif (dna(i).MeanIntensity*NumLab(i).Area < (g2-0.8*s2))
                        Red(NumLab(i).PixelIdxList) = 0;
                        Green(NumLab(i).PixelIdxList) = 1;
                        Blue(NumLab(i).PixelIdxList) = 0;
                        
                    
                    elseif (dna(i).MeanIntensity*NumLab(i).Area < (g2+3*s2))
                            Red(NumLab(i).PixelIdxList) = 0;
                            Green(NumLab(i).PixelIdxList) = 0;
                            Blue(NumLab(i).PixelIdxList) = 1;
                            
                            
                        else
                        
end
                        
        
        
end

mask = cat(3, Red, Green, Blue);
% figure, imshow(mask, []);

%%
cd ..
imwrite(mask, 'RGB.tif', 'compression','none'); 
imwrite(DAPI, 'DAPI.tif', 'compression','none'); 
imwrite(B, 'EU.tif', 'compression','none');

%%

