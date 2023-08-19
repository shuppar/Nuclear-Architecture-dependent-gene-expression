function [TotNucleoInt, AvgNucleoInt, AvgNucleoAre, TotNucPlsInt, AvgNuclPlsInt, numNucleoli] = NucleoliStat(QuantImage, Image, nuclearMask)



%%

A = Image; F = QuantImage;
A = imgaussfilt(A, 1);
A = imtophat(A, strel('Disk',29));
A = anisodiff(A, 5, 5, 0.25, 1);

%%

se = strel('disk', 5);
Io = imopen(A, se);
Ie = imerode(A, se);
Iobr = imreconstruct(Ie, A);
Ioc= imclose(Io, se);
Iobrd = imdilate(Iobr, se); Iobrcbr = imreconstruct(imcomplement(Iobrd), imcomplement(Iobr));
Iobrcbr = imcomplement(Iobrcbr);

%%

conli = stretchlim(Iobrcbr);
B = imadjust(Iobrcbr, conli, [0 1], 4); 
C = im2bw(B, graythresh(B));
C = bwareaopen(C, 160);

% figure, imshow(C); rect = getrect; C = imcrop(C, rect); F = imcrop(Image, rect);

%%

D = -bwdist(~C);
mask = imextendedmin(D, 3);   % The number (here 1) is very important... Decides the size of minima.
D = imimposemin(D, mask);

%%

E = watershed(D); 
L = C; L(E == 0) = 0; % figure, imshow(label2rgb(bwlabel(E),'jet',[.5 .5 .5]), []);
L = bwareaopen(L, 160);

%%

[L, numNucleoli] = bwlabel(L);
Stats = regionprops(L, F, 'Area', 'MeanIntensity');
Area = [Stats.Area]; MeanIntensity = [Stats.MeanIntensity];
Intensity = Area.*MeanIntensity;

if numNucleoli == 0
    
    AvgNucleoAre = 0;
    TotNucleoAre = 0;
    TotNucleoInt = 0;
    AvgNucleoInt = 0;
    
else
    
    AvgNucleoAre = mean(Area);
    TotNucleoAre = sum(Area);
    TotNucleoInt = sum(Intensity);
    AvgNucleoInt = TotNucleoInt/TotNucleoAre;
    
end
%% Nucleoplasm Intensity

[M, ~] = bwlabel(nuclearMask);
EUInt = regionprops(M, F, 'Area', 'MeanIntensity');

TotEUInt = [EUInt.Area].*[EUInt.MeanIntensity];
TotNucPlsInt = TotEUInt - TotNucleoInt;
AvgNuclPlsInt = TotNucPlsInt/((EUInt.Area)-TotNucleoAre);





