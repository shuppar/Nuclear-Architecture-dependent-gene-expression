function outims = LOG_filter(ims, hs, sigma)

% This generates the LOG filter itself.
% The bandwidth (here 1.5) may need to be changed depending
% on the pixel size of your camera and your microscope's
% optical characteristics.
if nargin < 1
    
        fprintf('\nError! Please give in a stack images, kernel size and sigma in the same order.\n');
    
    elseif nargin < 2
        hs = 7;     % these hs and sigma values are good for unbinned 60X image on OlyIX83
        sigma = 1;
    elseif nargin < 3
        fprintf('\nPlease give in a stack images, kernel size and sigma in the same order.\n');
    elseif nargin < 4
        hs = hs;
        sigma = sigma;
    else    
        fprintf('\nError! Please give in a stack images, kernel size and sigma in the same order.\n');    
        
end

H = -fspecial('log', hs, sigma);  % Changed to the present from -fspecial('log', 15, 1.5) on April1, 2017.

% The above change is because the RNA spots in our case have spread of <7 pixels diametrically.
% 1o = 2*round(3*1.5); 
% Instead of 15 it should be 10; this is beacuse 3*\sigma covers almost 98% of the Gaussian.

% Here, we amplify the signal by making the filter "3-D"
H = 1/3*cat(3,H,H,H);

% Apply the filter
outims = imfilter(ims,H,'replicate');

% Set all negative values to zero
outims(find(outims<0)) = 0;