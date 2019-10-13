%% ***************** Cell segmentation, to start with. ********************
% ***************A part of it is from Pedro Kostelec's blog: **************************
% ****blog.pedro.si/2014/04/basic-cell-segmentation-in-matlab.html?view=classic

%% A few things to note:

% So you might need to chage the disk size in strel in the second section.
% Calculate the minimum and maximum nuclear sizes and make changes
% accordingly in the section devoted to it.


%% Shuppar Script for Intensity plots.

clear all;
clc;
A = imread('DAPI1.tif');
DAPI = imread('DAPI1.tif');
% A = rgb2gray(A);
% A = A(:,:,1);
A = adapthisteq(A);                                 % Some local adjustments, to be able to detect the dimmer cells.
A = imclearborder(A);                               % Eliminating the objects on the borders.
A = wiener2(A, [10,10]);                            % Removing pixels smaller than the given size.

% figure(1), imshow(A), title('original DAPI image');
%figure(2), imshow(B), title('Original H2AX image');
%figure(3), imshow(C), title('Original Ki67 image');


%% Removing the problem of oversegmentation.

se = strel('disk', 20);
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
bw2 = imopen(bw1, strel('disk',20));                 % Morphological opening.
bw3 = bwareaopen(bw2, 100);                         % Removing cells with less than 100 pixels.

% figure(7), imshow(bw3), title('Image after a some processing');

bw3_perim = bwperim(bw3);   
overlay = imoverlay(A, bw3_perim, [1 .3 .3]);       % user-defined function by Stephen L Eddins.

% figure(8), imshow(overlay), title('Perimeter overlaid on the original');

%% Discovering Putative nucleus centroid.

maxs = imextendedmax(A, 5);
maxs = imclose(maxs, strel('disk',20));
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
        if chhotamota(i).Area <= 9000 ...                                                   % put in your own criteria for the area... read the starting note.
           && chhotamota(i).Area >= 1000 ...    
           && (chhotamota(i).Perimeter)^2/(4*pi*chhotamota(i).Area) <= 1.35 ...             % Circularity condition
           && (chhotamota(i).Perimeter)^2/(4*pi*chhotamota(i).Area) >= 0.65 ...
           && (4*chhotamota(i).Area)/(pi*(chhotamota(i).MajorAxisLength)^2) <= 1.35 ...     % Roundness.          
           && (4*chhotamota(i).Area)/(pi*(chhotamota(i).MajorAxisLength)^2) >= 0.65;    
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


%Cell counting:
[L, num] = bwlabel(L);

% num = numel(cm1);
fprintf('\n\nThe number of cells detected is: %d\n\n', num);
NumLab = regionprops(L, 'Centroid');


%% Overlaying the cells detected over the original image.

mask = im2bw(L, 1);
% figure(16), imshow(mask);
overlay2 = imoverlay(A, mask, [.5 .8 .3]);

% figure(17), imshow(overlay2), title('Mask');
% 
% %% Displaying the number over the detected nuclei
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

%% What does PixelIdxList do?
% B = A;
% for i = 1:(num/2)
%    idx = s(i).PixelIdxList;
%    B(idx) = 0;
%     
% end
% 
% figure(18), imshow(B);

%% Now try finding the intensity.


DI = zeros(num,1);      % DAPI intensity

for j = 1:num
    
    idx1 = d(j).MeanIntensity;
    idx2 = d(j).Area;
    DI(j) = idx1*idx2;
    f = fopen('Intensity.dat','a');  
    % fprintf(f,'%f\t%f\n', DI(j),d(j).MeanIntensity);
    fprintf(f,'%f\n', DI(j));

end

% filename = 'Intensity.dat';
% save(filename, 'I', '-ascii');

% hist(DI,10000);
 
 %% The end!
 
