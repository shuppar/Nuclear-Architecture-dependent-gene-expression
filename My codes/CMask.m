
function L = CMask(CellImage, NucleusImage, InterestNuc,varargin)

%Example: CellMask = CMask(cropPhal, nucmark, nucInt, 'diskR', 11, ...
% 'thresh','local', 'meth', 'm2', 'g1', 4, 'g2', 1.5);
% I have observed that g1 should be > 3.
% and g2 should be < 3.

%%
p = inputParser;

if nargin < 3
    error('MATLAB:narginchk:notEnoughInputs', 'Not enough input arguments.');
else 
    addParameter(p, 'thresh', 'global', ...
        @(t) (ischar(t) && ismember(t, {'global', 'local'})));
    addParameter(p, 's1', 'Mean', ...
        @(t) (ischar(t) && ismember(t, {'Mean', 'Gaussian'})));
    addParameter(p, 's2', 'Mean', ...
        @(t) (ischar(t) && ismember(t, {'Mean', 'Gaussian'})));
    addParameter(p, 'meth', 'm1', ...
        @(t) (ischar(t) && ismember(t, {'m1', 'm2'})));
    addParameter(p, 'diskR', 7, @(n) isinteger(uint8(n)));    % disk radius for imerode etc.
    addParameter(p, 'g1', 3, @(n) (isnumeric(n) && (n>=0)));  % gamma for global thresh
    addParameter(p, 'g2', 1.5,@(n) (isnumeric(n) && (n>=0))); % gamma for watershed
end
p.KeepUnmatched = true;

parse(p, varargin{:});
param = p.Results;
thresh = param.thresh;
Rad = param.diskR;
ga1 = param.g1;
ga2 = param.g2;
s1 = param.s1;
s2 = param.s2;
method = param.meth;

%% Initialization

% Works best with the images cropped with the cell in the centre...
% See boxlength function for more.

B = CellImage; DAPI = NucleusImage;
DAPI = uint16(anisodiff(DAPI, 17, 10, 0.25, 1));              % anisodiff(im, niter, kappa, lambda, option)

% This affects m2 quite a bit... play around with this if you will.
B = imgaussfilt(B, 0.35);
B = uint16(anisodiff(B, 17, 25, 0.25, 2));

%% Arguments for anisodiff:
%         im     - input image
%         niter  - number of iterations.
%         kappa  - conduction coefficient 20-100 ?
%         lambda - max value of .25 for stability
%         option - 1 Perona Malik diffusion equation No 1
%                  2 Perona Malik diffusion equation No 2
%
% Returns:
%         diff   - diffused image.
%
% kappa controls conduction as a function of gradient.  If kappa is low
% small intensity gradients are able to block conduction and hence diffusion
% across step edges.  A large value reduces the influence of intensity
% gradients on conduction.
%%
sz = size(B);
In = InterestNuc;
Blank = zeros(size(DAPI));

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
Atemp = imadjust(Aobrcbr, conli, [0 1], 0.5); 
C = im2bw(Atemp, graythresh(Atemp));

Amin = -bwdist(~C);
% bgm = watershed(Amin);
mask = imextendedmin(Amin, 3);   % The number (here 1) is very important... Decides the size of minima.
Amin = imimposemin(Amin, mask);
somet = watershed(Amin);
C(somet == 0) = 0;
% somet = imerode(C, strel('disk', 10));
somet = imdilate(C, strel('disk', 5));


%% Processing the cytosol image

se = strel('disk', Rad);
De = imerode(Dobrcbr, se);
Dobr = imreconstruct(De, Dobrcbr);
Dobrd = imdilate(Dobr, se);
Dobrcbr = imreconstruct(imcomplement(Dobrd), imcomplement(Dobr));
Dobrcbr = imcomplement(Dobrcbr);
Dobrcbr = uint16(anisodiff(Dobrcbr, 3, 5, 0.25, 1));

%% Background

        celi = stretchlim(imcomplement(Dobrcbr));
        Cetemp = imadjust(imcomplement(Dobrcbr), celi, [0 1], ga1);
        Imp = imcomplement(im2bw(Cetemp, graythresh(Cetemp)));
        Imp(somet==1) = 1;
        M = ~(bwareaopen(~Imp, 500));
        
        
        bgm = Blank;
        bgm1 = ~M; [t1, t2] = bwlabel(bgm1);
        bgm2 = regionprops(t1, 'Centroid');
        for i = 1:t2
        x = uint64(bgm2(i).Centroid); y = x(1); x = x(2);
        bgm(x, y) = 1;
        end
        clear bgm1; clear bgm2;
        bgm1 = imdilate(bgm, strel('disk', 5));
        bgm2 = imextendedmin(-bwdist(M), 1);
       
        dist = bwdist(M);
        distl = watershed(dist);
%         distl = watershed(imimposemin(uint8(dist), fgm));  % Not at all a good idea, this one!
        bgm = distl == 0;
        bgm(bgm1 == 1) = 1;
        bgm(bgm2 == 1) = 1; 
        clear bgm2; clear bgm1;


%% Adaptive threshold is must for Watershed lines.

switch method
    case 'm1'
       
        switch thresh
            case 'global'
                celi = stretchlim(imcomplement(Dobrcbr));
                Cetemp = imadjust(imcomplement(Dobrcbr), celi, [0 1], ga1);
                Imp = imcomplement(im2bw(Cetemp, graythresh(Cetemp)));
                Imp(somet==1) = 1;
                L = ~(bwareaopen(~Imp, 500));
            case 'local'
                celi = stretchlim((Dobrcbr));
                Dtemp = imadjust((Dobrcbr), celi, [0 1], ga1); 
                % L = im2bw(Dtemp, graythresh(Dtemp));
                threshp = adaptthresh(Dtemp, 0.7, 'Statistic', s1);  % you can try with Gaussian/Mean instead of Mean here. :)
                L = (imbinarize(Dtemp, threshp));
                L(somet==1) = 1;
                L = ~(bwareaopen(~L, 2000));
            otherwise
                fprintf('\n\nTheshold should be either "global" or "local"\n\n');
        end
        
        conli = stretchlim(Dobrcbr);
        Dtemp = imadjust(Dobrcbr, conli, [0 1], ga2); 
        % E = im2bw(Dtemp, graythresh(Dtemp));
        thresh2 = adaptthresh(Dtemp, 0.7, 'Statistic', s2);  % you can try with Gaussian/Mean instead of Mean here. :)
        E = imbinarize(Dtemp, thresh2);

        E(somet==1) = 1;
        E = ~(bwareaopen(~E, 2000));
        
        Dmin = -bwdist(~E);
        mask1 = imextendedmin(Dmin, 11, 8);   % The number (here 11) is very important... Decides the size of minima.

        % For Dmin1 below instead of C (which is nuclear mask) you can try using mask1, mask or somet.
        Dmin1 = imimposemin(Dmin, C|bgm);
        F = watershed(Dmin1);

        L(F == 0) = 0;
        
    case 'm2'

        hy = fspecial('sobel');
        hx = hy';
        Iy = imfilter(double(B), hy, 'replicate');
        Ix = imfilter(double(B), hx, 'replicate');
        gradmag = sqrt(Ix.^2 + Iy.^2);
        L = M;
        
        
% Some relic tries

%         fgm = Blank;
%         fgm1 = C; [t1, t2] = bwlabel(fgm1);
%         fgm2 = regionprops(t1, 'Centroid');
%         for i = 1:t2
%         x = uint64(fgm2(i).Centroid); y = x(1); x = x(2);
%         fgm(x, y) = 1;
%         end
%         clear fgm1; clear fgm2;
%         fgm = imdilate(fgm, strel('disk', 11));

        fgm = imerode(C, strel('disk', 7));

        gradmagp = imimposemin(gradmag, bgm|fgm);
        F = watershed(gradmagp);
        L(F==0) = 0;
        
    otherwise
end

%% Watershed for the final segmentation

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
       Lp(cellId(chun).PixelIdxList) = 1;
       Lp = imfill(Lp, 'holes');
       Lp = imopen(Lp, strel('disk', 2));
       Lp = bwareaopen(Lp, 300);
       L = Lp;
    else        
    end    
    
end

[~, cellnum] = bwlabel(L);
    
if cellnum == 1
    
    else
        fprintf('\n\nCould not find the cell of interest, sorry!\n\n');
end



end

%%

% figure, imshow(Dobrcbr, []); figure, imshow(label2rgb(bwlabel(L),'jet',[.5 .5 .5]), []);
% text(204, 204, '*','HorizontalAlignment', 'center','VerticalAlignment', 'middle', 'Color', 'r');

%% Display (For now commented) 

% test = label2rgb(bwlabel(L));
% figure, imshow(test);
% hold on;
% h = imshow(Dobrcbr, []);
% set(h, 'AlphaData', 0.5); hold off;

%% Some old failed tries for nostalgia

% hy = fspecial('sobel');
% hx = hy';
% Iy = imfilter(double(imclearborder(B)), hy, 'replicate');
% Ix = imfilter(double(imclearborder(B)), hx, 'replicate');
% gradmag = sqrt(Ix.^2 + Iy.^2);
% gradmag = uint16(anisodiff(gradmag, 17, 1, 0.25, 1));


% bw = im2bw(Dobrcbr, graythresh(Dobrcbr));
% Dist = bwdist(bw);
% DL = watershed(Dist);
% bgm = DL == 0;
% fgm = somet;
% gradmag2 = imimposemin(gradmag, fgm);
% M = watershed(gradmag2);
