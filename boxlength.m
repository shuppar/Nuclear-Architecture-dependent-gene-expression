
function bl = boxlength(labeledImage, objectNumber, k)


%% To find the box length of a square box around the nucleus of interest.

% This fucntion is used with CMask function to crop images for better cell segmentation.
% Boxlength is defined as the largest of the smallest (k) distances as measured from 
% the centroid of the nucleus of interest to the centroids of all other
% nuclei (sort of Hausdorff distance?). :)

n = objectNumber;
I = labeledImage;

[~, num] = bwlabel(I);
a = regionprops(I, 'Centroid');
a = cat(1, a.Centroid);
dist1 = zeros(num, 1);
bl2 = zeros(k,1);

for i = 1:num
    dist1(i) = sqrt((abs(a(n,1)-a(i,1)))*(abs(a(n,1)-a(i,1))) ...
        + (abs(a(n,2)-a(i,2)))*(abs(a(n,2)-a(i,2))));
    
end

dist = nonzeros(dist1);
bl1 = sortrows(dist, 1);

for i = 1:k
bl2(i) = bl1(i);
end

bl = max(bl2);

end