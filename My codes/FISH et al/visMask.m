function visMask = visMask(NuclearMask)


%%
alphaImage = im2double(NuclearMask);
NumLab = regionprops(NuclearMask, 'PixelIdxList');
[~, num] = bwlabel(alphaImage);
greenmask = alphaImage;
redmask = alphaImage;
bluemask = alphaImage;

for i = 1:num
    redmask(NumLab(i).PixelIdxList) = 0;
    greenmask(NumLab(i).PixelIdxList) = 255;
    bluemask(NumLab(i).PixelIdxList) = 0;
    
end


visMask = cat(3, redmask, greenmask, bluemask);

