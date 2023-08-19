%% function for 3d Nuclei Segmentation

function [NucSeg3d, Volume, Height, flist] = NMask3d(SNS, SNM, cutoff, sigma)

% SNS =  SignleNucStack... Stack of a single nucleus DAPI channel
% SNM = SingleNucMask2d... 2D Mask with segmented nucleus in the centre
% To get 2d mask, use NMask and then crop the rectangle containing the ...
% nucleus of interest
% cutoff is a 2d vector with values less than 1. This is to select only ...
% those images in the stack which have texture values more than the ...
% average of the complete stack multiplied by the cutoff... See the code ..
% for more info

if nargin < 3
    
        fprintf('\nPlease give in a stack of a single nucleus\n');
    
    elseif nargin < 4
        
        cA = SNS;
        cL = SNM;
        sz = 5;
        v1 = cutoff(1); % This is for initial images in the stack close to the coverslip.
        v2 = cutoff(2); % This is for images in the stack away from the coverslip
        
    elseif nargin < 5
         
        cA = SNS;
        sz = sigma;
        cL = SNM;
        v1 = cutoff(1); % This is for initial images in the stack close to the coverslip.
        v2 = cutoff(2); % This is for images in the stack away from the coverslip
        
    else
        
        fprintf('\nInvalid input.\n');
        
end

%%  Some initial processing for texture identification

[xc, yc, zc] = size(cA);
cL = imclearborder(cL);
% t1 = imgaussfilt3(cA, 0.5);
t2 = LOG_filter(cA, ceil(6*sz), sz);
t4 = t2; t5 = zeros(zc, 1);

%% Getting texture values

for i = 1:zc
    t4(:,:,i) = t2(:,:,i).*uint16(cL);
    temp = regionprops(cL, t4(:,:,i), 'MeanIntensity');
    t5(i) = temp(1).MeanIntensity; clear temp;
end
sel1 = v1*mean(t5); % Doing it (3*x)/4 is costlier than 0.75*x
sel2 = v2*mean(t5);

%% Removing out of focus images  from the stack.

% This is to select only those images in stack with nucleus present by ...
% removing out of focus planes.
ind = zeros(zc,1);
for i = 1:zc
    if t5(i) < sel1 && i < zc/2
    t4(:,:,i) = zeros(xc, yc);
    elseif t5(i) < sel2
       t4(:,:,i) = zeros(xc, yc); 
    else
        ind(i) = i;
    end
end

ind = ind(ind~=0); zi = length(ind);
Height = zi;

%% Finally making solid mask from the above texture images...

t3 = t4;
for i = 1:zc
t3(:,:,i) = imbinarize(t4(:,:,i)); % t3(:,:,i) = imclearborder(t3(:,:,i));
    for j = 1:2:4
        t3(:,:,i) = imdilate(t3(:,:,i), strel('disk', ceil(5*j))); % Have to set this 5 to a suitable number.
        t3(:,:,i) = imfill(t3(:,:,i), 'holes');
        t3(:,:,i) = imerode(t3(:,:,i), strel('disk', ceil(5*j)));
        t3(:,:,i) = bwareaopen(t3(:,:,i), 1500);
    end

end

NucSeg3d = logical(t3);

xy = cell(zi,1); z = xy;
Pixel_List = struct('xy', xy, 'z', z);


%% Volume and pixelList

Volume = regionprops(t3, 'Area');
Volume = Volume.Area;

for i = 1:zi
    
    if i == 1 || i == zi
        
        % Now remember PixelList gives values in the form ...
        % (x, y) which in matrix parlance correspond to (column, row) ...
        % That is x = column and y = row,.
        j = ind(i);
        tv1 = regionprops(t3(:,:,j), 'PixelList');
        Pixel_List(i).xy = tv1.PixelList; 
        Pixel_List(i).z = j; clear tv1;
        
    else
        j = ind(i);
        t5 = t3(:,:,j);
        t5 = t5 - imerode(t5, strel('disk', 1));
        tv1 = regionprops(t5, 'PixelList');
        Pixel_List(i).xy = tv1.PixelList;
        Pixel_List(i).z = j; clear tv1;
%         test1 = Pixel_List(5).xy;
%         test = zeros(size(t5)); test3 = sub2ind(size(t5), (test1(:,2)), (test1(:,1)));
%         test(test3) = 1; figure, imshow(test, []);
    end

end

%% Getting final list of pixels

len = 0;
for i = 1:zi
    len = len + length(Pixel_List(i).xy);
end

flist = zeros(len, 3);
dc = 0;

for i = 1:zi
    for j = 1:length(Pixel_List(i).xy)
        tv2 = Pixel_List(i).xy;
        flist(dc+j, 1) = tv2(j,1);
        flist(dc+j, 2) = tv2(j, 2);
        flist(dc+j,3) = ind(i);
    end
    dc = dc+j;
%     pause;
end

%%

end


