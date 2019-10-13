%% ***************** Cell segmentation, to start with. ******************************
% ***************A part of it is from Pedro Kostelec's blog: ************************
% ****blog.pedro.si/2014/04/basic-cell-segmentation-in-matlab.html?view=classic******

%% A few things to note (Important)...

% You will have to take a blank fluorescein image to correct for uneven illumination.
% So you might need to chage the disk size in strel in the second section.
% Calculate the minimum and maximum nuclear sizes and make changes
% accordingly in the section devoted to it.


%% Shuppar Script for Intensity plots.

clear all;
clc;
B = imread('ATRI1.tif');
% B = B(:,:,1);
C = imread('Hist1.tif');
% C = C(:,:,1);
A = imread('DAPI1.tif');
DAPI = imread('DAPI1.tif');

% A = rgb2gray(A);
% A = A(:,:,1);
%% Nonuniform illumination correction. Function saved in home/MatlabCode. Copy of Image Analyst's code. (August 2, 2016)

% A = BackgroundCorrect(A, BlankImage.tif);         % Take a blank fluorescein image.
% B = BackgroundCorrect(B, BlankImage.tif);
% C = BackgroundCorrect(C, BlankImage.tif);

%% Continues as before...

A = adapthisteq(A);                                 % Some local adjustments, to be able to detect the dimmer cells.
A = imclearborder(A);                               % Eliminating the objects on the borders.
A = wiener2(A, [10,10]);                            % Removing pixels smaller than the given size.

% figure(1), imshow(A), title('original DAPI image');
%figure(2), imshow(B), title('Original H2AX image');
%figure(3), imshow(C), title('Original Ki67 image');


%% Removing the problem of oversegmentation.

se = strel('disk', 18);
Ao = imopen(A, se);

% figure(4), imshow(Ao), title('opening Ao');

Ae = imerode(A, se);
Aobr = imreconstruct(Ae, A);
Aoc = imclose(Ao, se);
Aobrd = imdilate(Aobr, se);
Aobrcbr = imreconstruct(imcomplement(Aobrd), imcomplement(Aobr));
Aobrcbr = imcomplement(Aobrcbr);

% figure(5), imshow(Aobrcbr), title('Opening-closing by reconstruction');

A = Aobrcbr;



%% From here it continues as before.

bw = im2bw(A, graythresh(A));                                % graythresh was giving black image. Hence put a random value, works nevertheless for this case.

%figure(6), imshow(bw);

bw1 = imfill(bw, 'holes');
bw2 = imopen(bw1, strel('disk',18));                 % Morphological opening.
bw3 = bwareaopen(bw2, 100);                         % Removing cells with less than 100 pixels.

% figure(7), imshow(bw3), title('Image after a some processing');

bw3_perim = bwperim(bw3);   
overlay = imoverlay(A, bw3_perim, [1 .3 .3]);       % user-defined function by Stephen L Eddins.

% figure(8), imshow(overlay), title('Perimeter overlaid on the original');

%% Discovering Putative nucleus centroid.

maxs = imextendedmax(A, 5);
maxs = imclose(maxs, strel('disk',18));
maxs = imfill(maxs, 'holes');
maxs = bwareaopen(maxs, 100);
overlay1 = imoverlay(A, bw3_perim | maxs, [1 .3 .3]);

% figure(9), imshow(overlay1), title('Finding centroid');

%% Modifying the image so that the background pixels and the extended maxima pixles are forced to be the only minima in the image.

Jc = imcomplement(A);

% figure(10), imshow(Jc), title('Lets see');

A_mod = imimposemin(Jc, ~bw3 | maxs);

% figure(11), imshow(A_mod);

%% Watershed algorithm

L = watershed(A_mod);
labaledImage = label2rgb(L);
% figure(14), imshow(L);
% % L = L - bwareaopen(L, 6000);
% figure(15), imshow(labaledImage);
%% Counting cells and removing under- and over segmented nuclei.

[L, num1] = bwlabel(L);
chhotamota = regionprops(L, 'Area', 'PixelIdxList', 'Perimeter', 'MajorAxisLength');


cm1 = zeros(1, num1);
cm2 = zeros(1, num1);
s.pixels = zeros(1, num1);

%Detecting nuclei with optimum area.

for i = 1:num1
%         fprintf('%d\n', chhotamota(i).Area);
        if chhotamota(i).Area <= 35000 ...                                                  % put in your own criteria for the area... read the starting note.
           && chhotamota(i).Area >= 8000 ...    
           && (chhotamota(i).Perimeter)^2/(4*pi*chhotamota(i).Area) <= 1.25 ...             % Circularity condition
           && (chhotamota(i).Perimeter)^2/(4*pi*chhotamota(i).Area) >= 0.45 ...
           && (4*chhotamota(i).Area)/(pi*(chhotamota(i).MajorAxisLength)^2) <= 1.25 ...     % Roundness.          
           && (4*chhotamota(i).Area)/(pi*(chhotamota(i).MajorAxisLength)^2) >= 0.45;
           cm1(i) = i;
           
        else
            s(i).pixels = chhotamota(i).PixelIdxList;
            cm2(i) = i;
        end    
end

cm1 = cm1(cm1 ~= 0);
cm2 = cm2(cm2 ~= 0);


% % s(~cellfun('isempty',s))


% Removal of those pixels from our labaled image L which fall into nuclei
% having area not falling into our criteria.

num2 = numel(cm2);
for i = 1:num2
%     fprintf('%d\n', i);
    L(s(cm2(i)).pixels) = 0;
    
end

clear num1;
%Cell counting:
[L, num] = bwlabel(L);

% num = numel(cm1);
fprintf('\n\nThe number of cells detected is: %d\n\n', num);
NumLab = regionprops(L, 'Centroid', 'PixelIdxList');


%% Overlaying the cells detected over the original image.

mask = im2bw(L, 1);
% figure(16), imshow(mask);
overlay2 = imoverlay(A, mask, [.5 .8 .3]);

% figure(17), imshow(overlay2), title('Mask');
% 
% %% Displaying the numbers over the detected nuclei
% hold on
% for k = 1:num
%     numlab = NumLab(k).Centroid;
%     text(numlab(1), numlab(2), sprintf('%d', k), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
%     
% end
% hold off

%% Trying to get the list of pixels. (Not required here, though the functions listed are very useful.)

% S = regionprops('table', L, 'PixelList','PixelIdxList');  

%% Getting mean intensity of DAPI, H2AX and Ki67 for the detected cells.

d = regionprops(L, DAPI, 'Area','MeanIntensity');
h = regionprops(L, C, 'Area','MeanIntensity');
k = regionprops(L, B, 'Area','MeanIntensity');
hstd = regionprops(L, stdfilt(C), 'Area','MeanIntensity');
hent = regionprops(L, entropyfilt(C), 'Area','MeanIntensity');
kstd = regionprops(L, stdfilt(B), 'Area','MeanIntensity');
kent = regionprops(L, entropyfilt(B), 'Area','MeanIntensity');

%% Writing Data


DI = zeros(num,1);      % DAPI intensity
HI = zeros(num,1);      % H2AX intensity
KI = zeros(num,1);      % Ki67 Intensity

for j = 1:num
    
    idx1 = d(j).MeanIntensity;
    idx2 = d(j).Area;
    idx3 = h(j).MeanIntensity;
    idx4 = hent(j).MeanIntensity*h(j).Area;
    idx5 = hstd(j).MeanIntensity*h(j).Area;
    idx6 = k(j).MeanIntensity;
    idx7 = kent(j).MeanIntensity*k(j).Area;
    idx8 = kstd(j).MeanIntensity*k(j).Area;
    
    DI(j) = idx1*idx2;
    HI(j) = idx3*h(j).Area;
    KI(j) = idx6*k(j).Area;
    
    f = fopen('Intensity.dat','a');  
    fprintf(f,'%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n', DI(j), KI(j), idx6, idx7, idx8, HI(j), idx3, idx4, idx5);
    
      % DAPI_total, Ki67_total, Ki67_mean, Ki67_std_texture, Ki67_entropy_texture, H2A_Total, H2a_mean, H2a_std_texture, H2a_entropy_texture
    
end



%% Colouring Cells in different stages of cell cycle in different colours

if exist('G1_peak.dat', 'file')
 
peak1 = load('G1_peak.dat');    
peak1 = peak1(:,1);
Red = L;
Blue = L;
Green = L;

for i = 1:num
%         
        if (d(i).MeanIntensity*d(i).Area < 0.75*peak1)
            Red(NumLab(i).PixelIdxList) = 200;
            Green(NumLab(i).PixelIdxList) = 125;
            Blue(NumLab(i).PixelIdxList) = 0;
          
        else
             if (d(i).MeanIntensity*d(i).Area < 1.25*peak1)
                 Red(NumLab(i).PixelIdxList) = 255;
                 Green(NumLab(i).PixelIdxList) = 0;
                 Blue(NumLab(i).PixelIdxList) = 0;
                 
                 
             else
                    if (d(i).MeanIntensity*d(i).Area < 1.75*peak1)
                        Red(NumLab(i).PixelIdxList) = 0;
                        Green(NumLab(i).PixelIdxList) = 255;
                        Blue(NumLab(i).PixelIdxList) = 0;
                    
                    else
                        if (d(i).MeanIntensity*d(i).Area < 2.25*peak1)
                            Red(NumLab(i).PixelIdxList) = 0;
                            Green(NumLab(i).PixelIdxList) = 0;
                            Blue(NumLab(i).PixelIdxList) = 255;                        
                        else
                            
                            Red(NumLab(i).PixelIdxList) = 0;
                            Green(NumLab(i).PixelIdxList) = 125;
                            Blue(NumLab(i).PixelIdxList) = 200;
                        
                        end
                        
                    end
                    
              end
             
             
        end
        
        
end

Mask1 = cat(3, Red, Green, Blue);
clear Red; clear Blue, clear Green;
figure(19), imshow(C, []);
hold on
handle = imshow(Mask1);
hold off
set(handle, 'AlphaData', 0.3);

else

end





%% Putting texture measure on the detected nuclei

% figure(20), imshow(C, []);
 
% hold on
% for k = 1:num
%     numlab = NumLab(k).Centroid;
%     text(numlab(1), numlab(2), sprintf('%d', hent(k).MeanIntensity*h(k).Area ), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
%     
% end
% hold off



%%

% close all;

 
 %% The end!

 
 % loglog(c, b, '.', 'MarkerEdgeColor', [.5 0 .5], 'MarkerFaceColor', [.7 0 .7], 'MarkerSize', 5, 'Marker', 'o');       % For getting nice scatter plots
 % semilogy(a, b, '.', 'MarkerEdgeColor', [.5 .5 .0], 'MarkerFaceColor', [.7 .7 0], 'MarkerSize', 5, 'Marker', 'o');
 
