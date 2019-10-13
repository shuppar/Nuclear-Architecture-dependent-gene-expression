%% Shuppar Script for Nuclear Intensity Plots.

close all;
clear all;
clc;

fmain=sprintf('1');
fileFolder = fullfile('1');
dirDAPI = dir(fullfile(fileFolder, '*C0001.tif')); 
dirgene = dir(fullfile(fileFolder, '*C0003.tif'));
dirRNA = dir(fullfile(fileFolder, '*C0002.tif'));
[gene_S, ~] =  AvgNStack(fileFolder, dirgene);
[RNA_S, ~] =  AvgNStack(fileFolder, dirRNA);

[DAPI_S, ~] = AvgNStack(fileFolder, dirDAPI);

[xs, ys, zs] = size(DAPI_S); tot = xs*ys;

%% Nonuniform illumination correction. Function saved in home/MatlabCode. Copy of Image Analyst's code. (August 2, 2016)

BlankImage = imread('BlankImage.tif');         % Take a blank fluorescein image.
DAPI = uint16(mean(BackgroundCorrect(DAPI_S, BlankImage), 3));
gene = uint16(mean(BackgroundCorrect(gene_S, BlankImage),3));
Blank = zeros(size(DAPI)); [nrows, ncols] = size(DAPI);
RNA = uint16(max(BackgroundCorrect(RNA_S, BlankImage), [], 3));
% phal = phal + C/1.5;

% index = FokInd(DAPI_S); %To select best focused image from the stack

%% Rolling ball subtraction

% DAPI = imtophat(DAPI, strel('disk', 121));
% B = imtophat(B, strel('disk', 121));
% C = imtophat(C, strel('disk', 121));
% A = DAPI;

%% Subtracting Background

if exist('bg_Allele.dat', 'file')
    
    bagr = load('bg_Allele.dat');
    DB = bagr(1); GB = bagr(2); RB = bagr(3);
    
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
GB = max(gene_S, [], 3).*(M);
RB = max(RNA_S, [], 3).*(M);

DB = (mean(DB(:))*tot)/(Area); GB = (mean(GB(:))*tot)/(Area);
RB = (mean(RB(:))*tot)/(Area);
% CB = (mean(CB(:))*tot)/(Area);

f = fopen('bg.dat', 'w');  
fprintf(f,'%f\t%f\t%f\n', DB, GB, RB);


end

DAPI = DAPI - DB; D = gene - GB; RNA = RNA-RB;
DAPI(DAPI<1) = 0; D(D<1) = 0; RNA(RNA<1) = 0;
% C = C - CB; C(C<1) = 0;
A = DAPI;


%% Nuclear Mask

[L, ~] = NMask(A, 5, 2000, 200000, 0.25, 1.4, 0.25, 1.4);  % NMask(Image, DiskRadius, minArea, maxArea, minCircularity, maxCircularity, minRoundness, maxRoundness)
L = mimclearborder(L, 100); % Just for those cases where cells are to be segmented.
[~, num] = bwlabel(L);
% Circularity is bigger than Roundness (at least for an ellipse) but ...
% the way they are defined are always smaller than 1.
% Ellipse; because most nuclei are elliptic.

% num = numel(cm1);
fprintf('\n\nThe number of cells detected is: %d\n\n', num);
NumLab = regionprops(L, 'Centroid', 'PixelIdxList', 'BoundingBox', 'Area');
% visMask = visMask(L);
% names = {dirDAPI.name};
%  montage(reshape(DAPI_S,[xs, ys, 1, zs]),'DisplayRange',[]);

%% mRNA mask

% SpotCount(hFoc, 20, 121, 5, 100, 'hsize', 7, 'sigma', 1);
[Rspot, ~] = SpotCount(RNA_S, 10, 400, 10, 100, 'hsize', 9, 'sigma', 1.5);

%% Gene Position Enhancement

tic
disn = ceil(num*3);
% [gm, ~] = SpotCount(gene_S, 20, disn, 200, 1500, 'hsize', 24, 'sigma', 4); % Earlier thing
[gm, ~] = SpotCount(gene_S, 20, disn, 10, 800, 'hsize', 18, 'sigma', 3);
toc

%%
dna = regionprops(L, DAPI, 'Area', 'MeanIntensity', 'BoundingBox');
DI = zeros(num, 1); eRa = DI;
gnum = DI;
vol = DI; plist = cell(num,1);
gd = zeros(num,3);
gid = gd; gcd = gd;
gd2 = gd; gcd2 = gd;
clz = gd; % measures degree of closeness to the peripherry defined: gd/(gd+gcd).
clz2 = gd;
H_nuc = DI;
c_co = zeros(num,3);
exprsn = gd;

%% Writing Data
    

for j = 1:num
    idx1 = dna(j).MeanIntensity;
    idx2 = dna(j).Area;
    DI(j) = idx1*idx2;
    
%% Box length for cropping square for automated cell segmentation

%  Im1 = Im([x_start:x_end],[y_start:y_end],:); % For cropping image stack.

xc = ceil(dna(j).BoundingBox);
yc = shuppar(xc(2)-20); xw = dhuppar(xc(3)+50, ncols); 
yw = dhuppar(xc(4)+50, nrows); xc = (xc(1)-20);
xw = dhuppar(xc+xw, ncols);
yw = dhuppar(yc+yw, nrows);
nL1 = L(yc:yw, xc:xw, :); nL1 = imclearborder(nL1);
n1 = DAPI_S(yc:yw, xc:xw, :); [xr, yr, zr] = size(n1);
gc1 = gm(yc:yw, xc:xw, :);
[n3d, vol(j), H_nuc(j), plist{j}] = NMask3d(n1, nL1, [0.85, 0.85], 5);
vol(j) = vol(j)*0.07567*0.07567*0.3; % Volume of a voxel at 60x
eRa(j) = (0.239*vol(j))^(1/3);  % 3/(4*pi) = 0.239...To find radius of a sphere with the same volume.
H_nuc(j) = H_nuc(j)*0.3; % in micrometers
tatti = regionprops(n3d, 'Centroid');
c_co(j,:) = tatti.Centroid; clear tatti;
gc1 = gc1.*n3d;
[gcm,gnum(j)] = bwlabeln(gc1);
gcs = regionprops(gcm, 'Centroid', 'BoundingBox');
gc = zeros(3,3);

Rsp1 = Rspot(yc:yw, xc:xw, :);

%%

if gnum(j) == 3
    
%    figure, montage(reshape(n3d+gc1,[xr, yr, 1, zr]),'DisplayRange',[]);
%    pause;
    
    for i = 1:3
        gct = gcs(i).Centroid;
        gc(i,2) = gct(2); gc(i,3) = gct(3); gc(i,1) = gct(1);
        clear gct;
    end
    
    gy1 = gc(1,2); gz1 = gc(1,3); gx1 = gc(1,1);
    gy2 = gc(2,2); gz2 = gc(2,3); gx2 = gc(2,1);
    gy3 = gc(3,2); gz3 = gc(3,3); gx3 = gc(3,1);
    
    gid(j,1) = sqrt(((gx1-gx2)^2 + (gy1-gy2)^2)*0.005726 + 0.09*(gz1-gz2)^2);
    gid(j,2) = sqrt(((gx1-gx3)^2 + (gy1-gy3)^2)*0.005726 + 0.09*(gz1-gz3)^2);
    gid(j,3) = sqrt(((gx3-gx2)^2 + (gy3-gy2)^2)*0.005726 + 0.09*(gz3-gz2)^2);
    gid(j,:) = sort(gid(j,:));
 
    td = zeros(length(plist{j}), 1);
    
 %%
    for i = 1:3

        for k = 1:length(plist{j})
            td1 = plist(j); td1 = td1{1};
            td(k) = sqrt(((td1(k,1)-gc(i,1))^2 + (td1(k,2)-gc(i,2))^2)*0.005726 + 0.09*(td1(k,3)-gc(i,3))^2);
        end
        
        % Distance from the lamina in 3d
        gd(j,i) = min(td);
        
        % Distance from the Centroid in 3d
        gcd(j,i) = sqrt(((c_co(j,1)-gc(i,1))^2 + (c_co(j,2)-gc(i,2))^2)*0.005726 + 0.09*(c_co(j,3)-gc(i,3))^2);
        
        % Closeness metric in 3d
        clz(j,i) = gd(j,i)/(gd(j,i)+gcd(j,i));
        
        % Projecting 3d nuclei in 2d
        n2d = logical(max(n3d, [], 3));
        t2d1 = regionprops(n2d, 'Centroid');
        t2d1 = t2d1.Centroid;
        n2d1 = n2d - imerode(n2d, strel('disk', 1));
        t2d2 = regionprops(n2d1, 'PixelList');
        t2d2 = t2d2.PixelList;
        
        % Distance from the lamina in 2d
        gd2t = zeros(length(t2d2),1);
        for k = 1:length(t2d2)
            gd2t(k) = sqrt(((t2d2(k,1)-gc(i,1))^2 + (t2d2(k,2)-gc(i,2))^2)*0.005726);
        end
        gd2(j,i) = min(gd2t);
        % Distance from the centroid in 2d
        gcd2(j,i) = sqrt(((t2d1(1)-gc(i,1))^2 + (t2d1(2)-gc(i,2))^2)*0.005726);
        
        % Closeness metric in 2d
        clz2(j,i) = gd2(j,i)/(gd2(j,i)+gcd2(j,i));
        clear gd2t; clear t2d1; clear t2d2; clear n2d; clear n2d1;
        
        % Alleles Expressing? Yes == 1, No == 0.        
        [nrt, nct, ~] = size(Rsp1);
        xt = ceil(gcs(i).BoundingBox);
        rx = shuppar(xt(2)-1); ryw = dhuppar(xt(4)+2, nct); 
        rxw = dhuppar(xt(5)+2, nrt); ry = (xt(1)-1);
        rz = shuppar(xt(3)-1); rzw = dhuppar(xt(6)+2, zs);
        rxw = dhuppar(rx+rxw, nrt);
        ryw = dhuppar(ry+ryw, nct);
        rzw = dhuppar(rz+rzw, zs);
        Rt = Rsp1(rx:rxw, ry:ryw, rz:rzw);
        Rt = imclearborder(Rt);
%         figure, imshow(max(Rt, [], 3), []); pause;
        clear xt, clear rx, clear rxw, clear ryw, clear ry
        [~, nt] = bwlabeln(Rt); % clear Rt;
        if nt < 1
        else
            exprsn(j, i) = 1;
        end
        clear nt       
    end


%% relabeling indices so that they are in ascending order
    
     [clz(j,:), IT] = sort(clz(j,:));
     
     gcdt = gcd;
     gcdt(j,1) = gcd(j,IT(1)); gcdt(j,2) = gcd(j,IT(2)); gcdt(j,3) = gcd(j,IT(3));
     gcd = gcdt; clear gcdt
     
     gdt = gd;
     gdt(j,1) = gd(j,IT(1)); gdt(j,2) = gd(j,IT(2)); gdt(j,3) = gd(j,IT(3));
     gd = gdt; clear gdt
     
     gdt2 = gd2;
     gdt2(j,1) = gd2(j,IT(1)); gdt2(j,2) = gd2(j,IT(2)); gdt2(j,3) = gd2(j,IT(3));
     gd2 = gdt2; clear gdt2
     
     gcdt2 = gcd2;
     gcdt2(j,1) = gcd2(j,IT(1)); gcdt2(j,2) = gcd2(j,IT(2)); gcdt2(j,3) = gcd2(j,IT(3));
     gcd2 = gcdt2; clear gcdt2
     
     clzt2 = clz2;
     clzt2(j,1) = clz2(j,IT(1)); clzt2(j,2) = clz2(j,IT(2)); clzt2(j,3) = clz2(j,IT(3));
     clz2 = clzt2; clear clzt2
     
     exprsnt = exprsn;
     exprsnt(j,1) = exprsn(j,IT(1)); exprsnt(j,2) = exprsn(j,IT(2)); exprsnt(j,3) = exprsn(j,IT(3));
     exprsn = exprsnt; clear exprsnt
    
    clear idx1; clear idx2; clear idx3; clear idx4; clear avgInt; clear td;
    clear avgArea; clear numNuc; clear CytoRPL; clear CytoHPG; close all;
 
%%    
else
    
end

        f = fopen('Allele_exprsn.dat','a');  
        fprintf(f,'%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n', ...
        DI(j), vol(j), eRa(j), gnum(j), ...
        gd(j, 1), gd(j, 2), gd(j, 3), mean(gd(j, :)), ...
        gid(j, 1), gid(j, 2), gid(j, 3), mean(gid(j, :)), H_nuc(j), ...
        gcd(j,1), gcd(j,2), gcd(j,3), mean(gcd(j, :)), ...
        clz(j,1), clz(j,2), clz(j,3), mean(clz(j,:)), ...
        exprsn(j,1), exprsn(j,2), exprsn(j,3), ...
        gd2(j, 1), gd2(j, 2), gd2(j, 3), mean(gd2(j, :)), ...
        gcd2(j,1), gcd2(j,2), gcd2(j,3), mean(gcd2(j, :)), ...
        clz2(j,1), clz2(j,2), clz2(j,3), mean(clz2(j,:)));

end


%% Saving 3d stack images

test = gene_S(yc:yw, xc:xw, :);

for i = 1:zs
imwrite(n3d(:, :, i), 'segNuc_stack.tif', 'WriteMode', 'append',  'Compression','none');
imwrite(gc1(:, :, i), 'SegGene_stack.tif', 'WriteMode', 'append',  'Compression','none');
imwrite(test(:, :, i), 'gene_stack.tif', 'WriteMode', 'append',  'Compression','none');
imwrite(n1(:, :, i), 'nuc_stack.tif', 'WriteMode', 'append',  'Compression','none');

end



%%
