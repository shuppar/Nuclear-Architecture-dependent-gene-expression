function L = CMask1(CellImage, NucleusImage, InterestNuc)


%% Initialization

% Works best with the images cropped with the cell in the centre...
% See boxlength function for more.

B = CellImage; DAPI = NucleusImage;
DAPI = uint16(anisodiff(DAPI, 17, 1, 0.25, 1));              % anisodiff(im, niter, kappa, lambda, option)
B = uint16(anisodiff(B, 17, 1, 0.25, 1));
sz = size(B);
In = InterestNuc;
Bg = imcomplement(B);

%anisodiff is a function by Peter Kovesi from The University of Western
%Australia

%% Getting nuclear mask

Dobrcbr = B;
se = strel('disk', 4);
Ae = imerode(DAPI, se);
Aobr = imreconstruct(Ae, DAPI);
Aobrd = imdilate(Aobr, se);
Aobrcbr = imreconstruct(imcomplement(Aobrd), imcomplement(Aobr));
Aobrcbr = imcomplement(Aobrcbr);

% figure, imshow(B, []);
% rect = getrect;
% Aobrcbr = imcrop(Aobrcbr, rect); 
% Dobrcbr = imcrop(B, rect);   % This whole portion was just for the check.

conli = stretchlim(Aobrcbr);
Atemp = imadjust(Aobrcbr, conli, [0 1], 0.9); 
C = im2bw(Atemp, graythresh(Atemp));

Amin = -bwdist(~C);
mask = imextendedmin(Amin, 3);   % The number (here 1) is very important... Decides the size of minima.
Amin = imimposemin(Amin, mask);
somet = watershed(Amin);
C(somet == 0) = 0;
% somet = imerode(C, strel('disk', 10));
somet = imdilate(C, strel('disk', 5));


%% Processing the cytosol image

se = strel('disk', 5);
De = imerode(Dobrcbr, se);
Dobr = imreconstruct(De, Dobrcbr);
Dobrd = imdilate(Dobr, se);
Dobrcbr = imreconstruct(imcomplement(Dobrd), imcomplement(Dobr));
Dobrcbr = imcomplement(Dobrcbr);
Dobrcbr = uint16(anisodiff(Dobrcbr, 3, 5, 0.25, 1));



%% Creating binary cell mask from the image

% Bgobrcbr = B;
% se = strel('disk', 5);
% Bge = imerode(Bgobrcbr, se);
% Bgobr = imreconstruct(Bge, Bgobrcbr);
% Bgobrd = imdilate(Bgobr, se);
% Bgobrcbr = imreconstruct(imcomplement(Bgobrd), imcomplement(Bgobr));
% Bgobrcbr = imcomplement(Bgobrcbr);
% Bgobrcbr = uint16(anisodiff(Bgobrcbr, 3, 5, 0.25, 1));


celi = stretchlim(imcomplement(Dobrcbr));
Cetemp = imadjust(imcomplement(Dobrcbr), celi, [0 1], 5);
Imp = imcomplement(im2bw(Cetemp, graythresh(Cetemp)));

% Adaptive threshold not required for getting binary cell mask
% thresh1 = adaptthresh(Dtemp, 0.5, 'Statistic', 'Gaussian');
% Dtemp = imbinarize(Dtemp, thresh1);
% Imp = Dtemp;

Imp(somet==1) = 1;
Imp = ~(bwareaopen(~Imp, 500));
Imp = imopen(Imp, strel('disk', 5));

%% Adaptive threshold is must for Watershed lines.

conli = stretchlim(Dobrcbr);
Dtemp = imadjust(Dobrcbr, conli, [0 1], 3); 
% E = im2bw(Dtemp, graythresh(Dtemp));
thresh2 = adaptthresh(Dtemp, 0.7, 'Statistic', 'Mean');  % you can try with Gaussian/Mean instead of Mean here. :)
E = imbinarize(Dtemp, thresh2);

E(somet==1) = 1;
E = ~(bwareaopen(~E, 2000));

%% Watershed for the final segmentation

Dmin = -bwdist(~E);
mask1 = imextendedmin(Dmin, 11, 8);   % The number (here 11) is very important... Decides the size of minima.

% For Dmin1 below instead of C (which is nuclear mask) you can try using mask1, mask or somet.
Dmin1 = imimposemin(Dmin, C);
F = watershed(Dmin1);
L = Imp; L(F == 0) = 0;
L = ~(bwareaopen(~L, 500));
[L, ~] = bwlabel(L);
L = bwareaopen(L, 2000);

%% Getting the cell of interest using InterestNuc

[~, cellnum] = bwlabel(L);
cellId = regionprops(L, 'PixelIdxList');

nucInd = regionprops(In, 'Centroid');
nucInd = nucInd(1).Centroid; nucx = int32(nucInd(2)); nucy = int32(nucInd(1));
% Rows represent y-axis and columns x-axis... see for yourself.
centId = sub2ind(sz, nucx, nucy);
Lp = zeros(sz);

for i = 1:cellnum
   
    yes = ismember(centId, cellId(i).PixelIdxList);
    if yes == 1
       chun = i;
    else
    end    
    
end


Lp(cellId(chun).PixelIdxList) = 1;
L = Lp;

end

%%

% figure, imshow(Dobrcbr, []); figure, imshow(label2rgb(bwlabel(L),'jet',[.5 .5 .5]), []);
% text(1176, 100, sprintf('*'),'HorizontalAlignment', 'center','VerticalAlignment', 'middle', 'Color', 'r');

%% Display (For now commented) 

% test = label2rgb(bwlabel(L));
% figure, imshow(test);
% hold on;
% h = imshow(Dobrcbr, []);
% set(h, 'AlphaData', 0.5); hold off;

%% Some old failed tries for nostalgia

hy = fspecial('sobel');
hx = hy';
Iy = imfilter(double(imclearborder(B)), hy, 'replicate');
Ix = imfilter(double(imclearborder(B)), hx, 'replicate');
gradmag = sqrt(Ix.^2 + Iy.^2);
gradmag = uint16(anisodiff(gradmag, 17, 1, 0.25, 1));

% 
% bw = im2bw(Dobrcbr, graythresh(Dobrcbr));
% Dist = bwdist(bw);
% DL = watershed(Dist);
% bgm = DL == 0;
% fgm = somet;
Bgmp = imerode(imcomplement(Imp), strel('disk', 11));
Bgmp = (imcomplement(Imp));
fgm = C;
gradmag2 = imimposemin(gradmag, fgm);
M = watershed(gradmag2);

