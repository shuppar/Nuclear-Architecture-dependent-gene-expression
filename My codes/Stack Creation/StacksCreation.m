%% Max and Average stack


clear;
close all;
clc;
fmain=sprintf('1');
%fimage=sprintf('%s-001.tif_Files', fmain);
fileFolder = fullfile('1'); %fimage removed
dirDAPI = dir(fullfile(fileFolder, '*C0001.tif'));  % Protein
dirkisi = dir(fullfile(fileFolder, '*C0002.tif'));
dirhist = dir(fullfile(fileFolder, '*C0003.tif'));
%[A, Max] = MaxNStack(fileFolder, dirProt);
[bDAPI, AvgDAPI] = AvgNStack(fileFolder, dirDAPI);
[bkisi, Avgkisi] = AvgNStack(fileFolder, dirkisi);
[bhist, Avghist] = AvgNStack(fileFolder, dirhist);
%imwrite(Max, 'MaxStack/1.tif', 'tif');
imwrite(AvgDAPI, 'AvgStack/DAPI1.tif', 'tif');
imwrite(Avgkisi, 'AvgStack/kisi1.tif', 'tif');
imwrite(Avghist, 'AvgStack/hist1.tif', 'tif');
