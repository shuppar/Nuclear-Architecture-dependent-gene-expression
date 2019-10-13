
%% Modified imclearborder for clearing a thick border in an image.


function O = mimclearborder(binaryImage, borderSize)


%%
b = borderSize;
I = binaryImage;
[r, c] = size(I);
w = (c-1) - 2*b ; h = (r-1) - 2*b;
CI = imcrop(I, [b+1, b+1, w, h]);
I = imclearborder(CI);
O = logical(padarray(I, [b b],'both'));

end

% The below is to check the alignment of the input and output images.
% choosing the coordinates (the first two entries) one can see that the
% above function gives perfect alignment.
% text(1176, 100, sprintf('*'),'HorizontalAlignment', 'center','VerticalAlignment', 'middle', 'Color', 'r');