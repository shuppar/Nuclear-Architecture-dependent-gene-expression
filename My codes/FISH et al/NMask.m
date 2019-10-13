
function [L, num] = NMask(Image, DiskRadius, minArea, maxArea, minCircularity, maxCircularity, minRoundness, maxRoundness)



%% Some Image processing to deal with over segmentation


Q = anisodiff(Image, 7, 5, 0.25, 2);              % anisodiff(im, niter, kappa, lambda, option)
Q = anisodiff(Q, 7, 5, 0.25, 1);              % anisodiff(im, niter, kappa, lambda, option)

%anisodiff is a function by Peter Kovesi from The University if Western
%Australia

A = uint16(Q);
conli = stretchlim(A);
A = imadjust(A, conli, [0 1], 0.8);          % Correction factor works fine in [0.7, 0.9].
A = wiener2(A, [3*DiskRadius, 3*DiskRadius]);                          % Removing pixels smaller than the given size.
se = strel('disk', DiskRadius);
Ae = imerode(A, se);
Aobr = imreconstruct(Ae, A);
Aobrd = imdilate(Aobr, se);
Aobrcbr = imreconstruct(imcomplement(Aobrd), imcomplement(Aobr));
Aobrcbr = imcomplement(Aobrcbr);
bw = im2bw(Aobrcbr, graythresh(Aobrcbr));

%% Added 11 August 2018... Earlier version was just
% without this. The next section coule be uncommented and ...
% this section commented to go back to earlier version.

Amin = -bwdist(~bw);
mask = imextendedmin(Amin, DiskRadius/2);   % Decides the size of minima.
Amin = imimposemin(Amin, mask);
somet = watershed(Amin); % Watershed lines
L = bw;
L(somet == 0) = 0;
L = imclearborder(L, 8);
[~, num1] = bwlabel(L);

%% Now this will be nothing like Pedro's blog. % [Commented on August 11 2018.]

% hy = fspecial('sobel');
% hx = hy';
% Iy = imfilter(double(imclearborder(Image)), hy, 'replicate');
% Ix = imfilter(double(imclearborder(Image)), hx, 'replicate');
% gradmag = sqrt(Ix.^2 + Iy.^2);
% 
% 
% fgm = imregionalmax(Aobrcbr);
% fgm = imclose(fgm, strel('disk', 2));
% fgm = imerode(fgm, strel('disk', 2));
% fgm = bwareaopen(fgm, 3*DiskRadius);
% Dist = bwdist(bw);
% DL = watershed(Dist);
% bgm = DL == 0;
% gradmag2 = imimposemin(gradmag, bgm | fgm);
% L = watershed(gradmag2);
% 
% [L, num1] = bwlabel(L);
% L = imclearborder(L, 8);


%% Counting and removing under- and over-segmented nuclei.

chhotamota = regionprops(L, 'Area', 'PixelIdxList', 'Perimeter', 'MajorAxisLength');


for i = 1:num1
    
        if chhotamota(i).Area <= maxArea ...                                                          % put in your own criteria for the area... read the starting note.
           && chhotamota(i).Area >= minArea ...    
           && (4*pi*chhotamota(i).Area)/(chhotamota(i).Perimeter)^2 <= maxCircularity ...             % Circularity condition
           && (4*pi*chhotamota(i).Area)/(chhotamota(i).Perimeter)^2 >= minCircularity ...
           && (4*chhotamota(i).Area)/(pi*(chhotamota(i).MajorAxisLength)^2) <= maxRoundness ...       % Roundness.          
           && (4*chhotamota(i).Area)/(pi*(chhotamota(i).MajorAxisLength)^2) >= minRoundness;
           
        else
              L(chhotamota(i).PixelIdxList) = 0;
        end    
end

L = imopen(L, strel('disk', 7));       %April 7, 2017; to smooth out the edges.
[~, num] = bwlabel(L);




%% dhuppar function

function out = dhuppar(nrows, x)

if nrows - x > 0
    out = x;
else
    out = nrows;
    
end

%% shuppar function

function out = shuppar(x)

if x > 0
    out = x;
else
    out = 1;
    
end


%% Aniso Diff function


% ANISODIFF - Anisotropic diffusion.
%
% Usage:
%  diff = anisodiff(im, niter, kappa, lambda, option)
%
% Arguments:
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
%
% lambda controls speed of diffusion (you usually want it at a maximum of
% 0.25)
%
% Diffusion equation 1 favours high contrast edges over low contrast ones.
% Diffusion equation 2 favours wide regions over smaller ones.

% Reference: 
% P. Perona and J. Malik. 
% Scale-space and edge detection using ansotropic diffusion.
% IEEE Transactions on Pattern Analysis and Machine Intelligence, 
% 12(7):629-639, July 1990.
%
% Peter Kovesi  
% School of Computer Science & Software Engineering
% The University of Western Australia
% pk @ csse uwa edu au
% http://www.csse.uwa.edu.au
%
% June 2000  original version.       
% March 2002 corrected diffusion eqn No 2.

function diff = anisodiff(im, niter, kappa, lambda, option)

if ndims(im)==3
  error('Anisodiff only operates on 2D grey-scale images');
end

im = double(im);
[rows,cols] = size(im);
diff = im;
  
for i = 1:niter
%  fprintf('\rIteration %d',i);

  % Construct diffl which is the same as diff but
  % has an extra padding of zeros around it.
  diffl = zeros(rows+2, cols+2);
  diffl(2:rows+1, 2:cols+1) = diff;

  % North, South, East and West differences
  deltaN = diffl(1:rows,2:cols+1)   - diff;
  deltaS = diffl(3:rows+2,2:cols+1) - diff;
  deltaE = diffl(2:rows+1,3:cols+2) - diff;
  deltaW = diffl(2:rows+1,1:cols)   - diff;

  % Conduction

  if option == 1
    cN = exp(-(deltaN/kappa).^2);
    cS = exp(-(deltaS/kappa).^2);
    cE = exp(-(deltaE/kappa).^2);
    cW = exp(-(deltaW/kappa).^2);
  elseif option == 2
    cN = 1./(1 + (deltaN/kappa).^2);
    cS = 1./(1 + (deltaS/kappa).^2);
    cE = 1./(1 + (deltaE/kappa).^2);
    cW = 1./(1 + (deltaW/kappa).^2);
  end

  diff = diff + lambda*(cN.*deltaN + cS.*deltaS + cE.*deltaE + cW.*deltaW);

%  Uncomment the following to see a progression of images
%  subplot(ceil(sqrt(niter)),ceil(sqrt(niter)), i)
%  imagesc(diff), colormap(gray), axis image

end
%fprintf('\n');


%% Display (For now commented) 

% test = label2rgb(bwlabel(L));
% figure, imshow(test);
% hold on;
% h = imshow(Dobrcbr, []);
% set(h, 'AlphaData', 0.5); hold off;