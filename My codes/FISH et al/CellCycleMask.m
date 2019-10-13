function mask = CellCycleMask(DAPI, L, peak1)
%%

g1 = peak1(1,1);
s1 = peak1(1,3);
g2 = peak1(1,2);
s2 = peak1(1,4);
Red = zeros(size(L));
Blue = Red; Green = Red;
NumLab = regionprops(L, DAPI, 'PixelIdxList', 'MeanIntensity', 'Area');
[~, num] = bwlabel(L);

%%


for i = 1:num

        if (NumLab(i).MeanIntensity*NumLab(i).Area < (g1-3*s1))  
            Red(NumLab(i).PixelIdxList) = 255/255;
            Green(NumLab(i).PixelIdxList) = 255/255;
            Blue(NumLab(i).PixelIdxList) = 178/255;
            fprintf('G1-\n');

        else
             if (NumLab(i).MeanIntensity*NumLab(i).Area < (g1+2.0*s1))
                 Red(NumLab(i).PixelIdxList) = 254/255;
                 Green(NumLab(i).PixelIdxList) = 204/255;
                 Blue(NumLab(i).PixelIdxList) = 92/255;
                 fprintf('G1\n');
                 
             else
                    if (NumLab(i).MeanIntensity*NumLab(i).Area < (g2-0.8*s2))
                        Red(NumLab(i).PixelIdxList) = 253/255;
                        Green(NumLab(i).PixelIdxList) = 141/255;
                        Blue(NumLab(i).PixelIdxList) = 60/255;
                        fprintf('S\n');
                    
                    else
                        if (NumLab(i).MeanIntensity*NumLab(i).Area < (g2+3*s2))
                            Red(NumLab(i).PixelIdxList) = 240/255;
                            Green(NumLab(i).PixelIdxList) = 59/255;
                            Blue(NumLab(i).PixelIdxList) = 32/255;
                            fprintf('G2\n');
                            
                        else
                            
                                 Red(NumLab(i).PixelIdxList) = 189/255;
                                 Green(NumLab(i).PixelIdxList) = 0/255;
                                 Blue(NumLab(i).PixelIdxList) = 38/255;
                                 fprintf('G2+');
                        
                        end
                        
                    end
                    
              end
             
             
        end
        
        
end

mask = cat(3, Red, Green, Blue);
figure, imshow(mask, []);


%%




