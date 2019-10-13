 

function [L, num] = SpotCount(Image_Stack, SpotRadius, DiscardNum, minArea, maxArea, varargin)

% For example: [hFoc, ~] = SpotCount(hFoc, 20, 121, 5, 100, 'hsize', 7, 'sigma', 1);

%%

p = inputParser;

if nargin < 3
    error('MATLAB:narginchk:notEnoughInputs', 'Not enough input arguments.');
else 
    addParameter(p, 'hsize', 7, @(n) isinteger(uint8(n)));    % disk radius for imerode etc.
    addParameter(p, 'sigma', 1, @(n) (isnumeric(n) && (n>=0)));
end
p.KeepUnmatched = false;

parse(p, varargin{:});
param = p.Results;
hsz = param.hsize;
sig = param.sigma;


%%

Stack = Image_Stack;
sr = SpotRadius;
Stack = imtophat(Stack, strel('disk', sr));

%%
Stack = double(Stack);
LF = -fspecial('log', hsz, sig);
LF =  1/3*cat(3,LF,LF,LF);
Stack = imfilter(Stack, LF, 'replicate');
Stack(find(Stack<0)) = 0;
Stack = Stack/max(Stack(:));

temp = max(Stack, [], 3);
temp1 = temp;
[nrows, ncols] = size(temp);


%% 

dn = DiscardNum;
threshGv = zeros(dn, 1);


for i = 1:dn
   
    [~, idx] = max(temp1(:));
    [x, y]   = ind2sub(size(temp1), idx);
    gtbox    = temp1(shuppar((x-20)):dhuppar(nrows, (x+20)), shuppar((y-20)):dhuppar(ncols, (y+20)));
    temp2    = max(gtbox(:));
    gtbox    = gtbox/temp2;
    threshGv(i) = graythresh(gtbox)*temp2; % temp1 added July29_2017;
    temp1(shuppar((x-20)):dhuppar(nrows, (x+20)), shuppar((y-20)):dhuppar(ncols, (y+20))) = 0;
    
end

thresh = mean(threshGv);
L   = Stack > thresh;
% figure, imshow(L(:,:,8), []);
% thresh

%% Removing false spots

L   = bwareaopen(L, minArea);

[~, num1] = bwlabeln(L);
temp2 = regionprops(L, 'Area', 'PixelIdxList');

for i = 1:num1
    if temp2(i).Area <= maxArea
       
    else
          L(temp2(i).PixelIdxList) = 0;
    end
end

[~, num] = bwlabeln(L);
 


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


