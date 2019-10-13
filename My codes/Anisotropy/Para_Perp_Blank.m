%% Shuppar Script Trial

I = imread('Focused/Image_1904.tif'); 
% figure, imshow(I, []);
P = imread('Image_Perp.tif');      % Mask images generated in ImageJ... Next time don't forget to get blank image.
L = imread('Image_PAra.tif');
P = (P>30);
L = (L>30);
for i = 1:3;
P = imopen(P, strel('square', 2));
L = imopen(L, strel('square', 2));
P = imfill(P, 'holes');
L = imfill(L, 'holes');
end
EP = edge(P,'canny');
EPd = imdilate(EP,strel('disk',10));
Pfilt = imfilter(P,fspecial('gaussian'));
P(EPd) = Pfilt(EPd);
LP = edge(L,'canny');
ELd = imdilate(EL,strel('disk',10));
Lfilt = imfilter(L,fspecial('gaussian'));
L(ELd) = Lfilt(ELd);
figure, imshow(P);
figure, imshow(L);

%% Some hints
% %Standard IPT Image
% I = imread('cameraman.tif');
% %Its edges
% E = edge(I,'canny');
% %Dilate the edges
% Ed = imdilate(E,strel('disk',2));
% %Filtered image
% Ifilt = imfilter(I,fspecial('gaussian'));
% %Use Ed as logical index into I to and replace with Ifilt
% I(Ed) = Ifilt(Ed);

%% hint 2
