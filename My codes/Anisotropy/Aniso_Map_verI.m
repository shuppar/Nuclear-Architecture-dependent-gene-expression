%% Shuppar Script for Anisotropy Map. %%a

close all;
clear all;
clc;

%% Mask and crop.

I = imread('Image_1527.tif'); 
% figure(8), imshow(I, []);

%% Parallel and Perpendicular Channel Mask
% 
P = imread('Image_Perp.tif');      
L = imread('Image_Para.tif');
P = (P>30);
L = (L>30);
for i = 1:3;
P = imopen(P, strel('square', 2));
L = imopen(L, strel('square', 2));
P = imfill(P, 'holes');
L = imfill(L, 'holes');
end
EP = edge(P,'canny');
EPd = imdilate(EP,strel('disk',10));
Pfilt = imfilter(P,fspecial('gaussian'));
P(EPd) = Pfilt(EPd);
EL = edge(L,'canny');
ELd = imdilate(EL,strel('disk',10));
Lfilt = imfilter(L,fspecial('gaussian'));
L(ELd) = Lfilt(ELd);
% figure(1), imshow(P);
% figure(2), imshow(L);

%% Generating Mask

A = bwboundaries(P);
B = bwboundaries(L);
MP = A{1};
ML = B{1};

xminp = min(MP(:,1)); yminp = min(MP(:,2));
xmaxp = max(MP(:,1)); ymaxp = max(MP(:,2));
widthp = (xmaxp - xminp); lengthp = (ymaxp - yminp);
C = imcrop(I, [yminp xminp lengthp widthp]);

% bounding box using regionprops

% SP = regionprops(P, 'BoundingBox');
% SL = regionprops(L, 'BoundingBox');
% C = imcrop(I, SP(1).BoundingBox);
% D = imcrop(I, SL(1).BoundingBox);


xminl = min(ML(:,1)); yminl = min(ML(:,2));
xmaxl = max(ML(:,1)); ymaxl = max(ML(:,2));
widthl = (max(ML(:,1)) - xminl); lengthl = (max(ML(:,2)) - yminl);
D = imcrop(I, [yminl xminl lengthp widthp]);


% figure(3), imshow(C, []);
% figure(4), imshow(D, []);

%% Aligning the cropped images.


tformEstimate = imregcorr(C, D, 'similarity');
Rfixed = imref2d(size(D));
C = imwarp(C, tformEstimate, 'OutputView', Rfixed);

E = imfuse(C,D,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);

figure(5), imshow(E);


%% Anisotropy formulae.  G Factor for 100x = 0.98.

C = im2double(C);
D = im2double(D);
Anisn = im2double(C - (0.98*D));
Anisd = im2double(C + (2*0.98*D));
Anis = im2double(imdivide(Anisn, Anisd));
colormap('hot');

figure(7), imagesc(Anis);
colorbar;
[m,n] = size(Anis);
figure(9), plot(Anis(floor((m/2)),:), '-');


%% Junk Lines

% figure, imshowpair(D, C, 'montage');
% figure(6), imshowpair(D, C, 'falsecolor');
