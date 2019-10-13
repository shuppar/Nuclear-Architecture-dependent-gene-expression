%% Segmentation Function
function [N, num] = CytoSeg(AvgDAPI, MaxCytoImage, NucleiDiskRadius, CellDiskRadius)



%% 
Z = anisodiff(AvgDAPI, 5, 20, 0.25, 2);
C = anisodiff(MaxCytoImage, 5, 20, 0.25, 2);


%% Nuclei seg


hy = fspecial('sobel');
hx = hy';
seN = strel('disk', NucleiDiskRadius);
A = uint16(Z);

%%

conli = stretchlim(A);
A = imadjust(A, conli, [0 1], 1.3);                 % April 1, 2017. Gamma Correction.

%%

A = imclearborder(A, 8);                            % Eliminating the objects on the borders.

% April 8, 2017: It seems the above one - omclearborder - is really an important step; removing
% which disrupts the whole segmentation with just nuclei remaining in the whole segmented cell.

A = wiener2(A, [50,50]);                            % Removing pixels smaller than the given size.
Ae = imerode(A, seN);
Aobr = imreconstruct(Ae, A);
Aobrd = imdilate(Aobr, seN);
Aobrcbr = imreconstruct(imcomplement(Aobrd), imcomplement(Aobr));
Aobrcbr = imcomplement(Aobrcbr);
Aobrcbr = im2double(Aobrcbr);
Aobrcbr = Aobrcbr/max(Aobrcbr(:));

IyD = imfilter(double(Z), hy, 'replicate');
IxD = imfilter(double(Z), hx, 'replicate');
gradmagD = sqrt(IxD.^2 + IyD.^2);

fgmD = imregionalmax(Aobrcbr);
fgmD = imclose(fgmD, strel('disk', 4));
fgmD = imerode(fgmD, strel('disk', 4));
fgmD = bwareaopen(fgmD, 100);
bwD = im2bw(Aobrcbr, graythresh(Aobrcbr));
DistD = bwdist(bwD);
DLD = watershed(DistD);
bgmD = DLD == 0;
gradmag2D = imimposemin(gradmagD, bgmD | fgmD);
L = watershed(gradmag2D);
[L,~] = bwlabel(L);
L(L==1) = 0;

%% Another something...

B = C;
B = imtophat(B, strel('disk', 256));  % you might have to play around with the disk radius.


%% Removing the problem of oversegmentation.

seC = strel('disk', CellDiskRadius);
Be = imerode(B, seC);
Bobr = imreconstruct(Be, B);

%%
Bobrd = imdilate(Bobr, seC);
Bobrcbr = imreconstruct(imcomplement(Bobrd), imcomplement(Bobr));
Bobrcbr = imcomplement(Bobrcbr);
Bobrcbrs = imgaussfilt(Bobrcbr, 1.5);
Bobrcbrs = im2double(Bobrcbrs);
Bobrcbrs = Bobrcbrs/max(Bobrcbrs(:));

%% Now this will be nothing like Pedro's blog.

Iy = imfilter(double(Bobrcbrs), hy, 'replicate');
Ix = imfilter(double(Bobrcbrs), hx, 'replicate');
gradmag = sqrt(Ix.^2 + Iy.^2);


%% Meyer Watershed, if you want to try anytime!

% fgm = imregionalmax(Dobrcbr);
% fgm = imclose(fgm, strel('disk', 4));
% fgm = imerode(fgm, strel('disk', 4));
% fgm = bwareaopen(fgm, 500);

% Dobrcbr = im2double(Dobrcbr);
% Dobrcbr = (Dobrcbr)/max(Dobrcbr(:));
% bw = im2bw(Dobrcbr, graythresh(Dobrcbr));
% Dist = bwdist(bw);
% DL = watershed(Dist);
% bgm = DL ==0;

% gradmag2 = imimposemin(gradmag, L|bgmD);   % Try replacing L with fgm.
% M = watershed(gradmag2);
% [M, num1] = bwlabel(M);

Cythresh = graythresh(Bobrcbrs);
N = IdentifySecPropagateSubfunction(L, im2double(B), imfill(im2bw(Bobrcbrs, 0.75*Cythresh), 'holes'), 0.05);  %         PropagatedImage = IdentifySecPropagateSubfunction(PrelimPrimaryLabelMatrixImage,OrigImage,ThresholdedOrigImage,RegularizationFactor)
N = imopen(N, strel('disk', 20));           % April 7, 2017; to smooth out edges.
num = max(N(:));







%% Nlaplace mask
% 
% H = ones(41);
% H(ceil((41^2)/2)) = 1 - 41^2;
% gradmag = imfilter(D, H, 'replicate');
% gradmag(find(gradmag<0)) = 0;



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



