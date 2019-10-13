%% Shuppar Script for Anisotropy Map. %%a

close all;
clear all;
clc;

%% Mask and crop.

I = imread('/home/shuppar/Downloads/Aniso/May 4, 2016 (Aniso_EGFP_and_H2B_transient)/EFGP/Image_2168.tif'); 
% figure(8), imshow(I, []);

%% Parallel and Perpendicular Channel Mask

P = imread('/home/shuppar/Downloads/Aniso/Image_Perp.tif');      
L = imread('/home/shuppar/Downloads/Aniso/Image_Para.tif');
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

%% Crop method I (Either use this or the one below)...  Remember C is Perpendicular and D is parallel

% xminp = min(MP(:,1)); yminp = min(MP(:,2));
% xmaxp = max(MP(:,1)); ymaxp = max(MP(:,2));
% widthp = (xmaxp - xminp); lengthp = (ymaxp - yminp);
% C = imcrop(I, [yminp xminp lengthp widthp]);    

% xminl = min(ML(:,1)); yminl = min(ML(:,2));translation
% xmaxl = max(ML(:,1)); ymaxl = max(ML(:,2));
% widthl = (max(ML(:,1)) - xminl); lengthl = (max(ML(:,2)) - yminl);
% D = imcrop(I, [yminl xminl lengthl widthl]);   



%% Croping method II (Either use this or the one above)... Remember C is Perpendicular and D is parallel

% bounding box using regionprops
SP = regionprops(P, 'BoundingBox');
SL = regionprops(L, 'BoundingBox');
C = imcrop(I, SP(1).BoundingBox);
D = imcrop(I, SL(1).BoundingBox);



%% Visualising the cropping.

% figure(3), imshow(C, []);    
% figure(4), imshow(D, []);

figure(10), imshow(I, []); hold on; visboundaries(A); hold on; visboundaries(B);

%% Aligning the cropped images.


tformEstimate = imregcorr(C, D, 'rigid');    % Try wuth different parameters: rigid, translation or similarity. 'Rigid' does not resize the images, while 'similarity' does.
Rfixed = imref2d(size(D));
C = imwarp(C, tformEstimate, 'OutputView', Rfixed);

E = imfuse(D,C,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);

figure(5), imshow(E);


%% Anisotropy formulae.  G Factor for 100x = 0.93.  ...  Remember C is Perpendicular and D is parallel

C = im2double(C);
D = im2double(D);
Anisn = im2double(D - (0.93*C));      
Anisd = im2double(D + (2*0.93*C));
Anis = im2double(imdivide(Anisn, Anisd));

figure(7), imagesc(Anis);
colormap jet;

cb = colorbar;
cb.Limits = [0 0.35];
[m,n] = size(Anis);
figure, plot(Anis(floor((m/2)),:), '-');


%% Junk Lines

% figure, imshowpair(D, C, 'montage');
% figure(6), imshowpair(D, C, 'falsecolor');

