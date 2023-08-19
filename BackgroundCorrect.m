%% Trying out clearing background for grey scale image from Image Analysts' generalised code. It's merely a copy of his code, with unwanted parts removed.


function correctedImage = BackgroundCorrect(inputImage, nonuniformBackgroundImage)
    
            [~, ~, zs] = size(inputImage);
            correctedImage = zeros(size(inputImage));
            inputImage = double(inputImage);
            nonuniformBackgroundImage = double(nonuniformBackgroundImage);
            
            
            for i = 1:zs
                noiselessImage = nonuniformBackgroundImage;
                maxValue = max(max(noiselessImage));
                modeledBackgroundImage = noiselessImage / maxValue;
                correctedImage(:,:,i) = (inputImage(:,:,i) ./ modeledBackgroundImage);
            end
            
            correctedImage = uint16(correctedImage);
end