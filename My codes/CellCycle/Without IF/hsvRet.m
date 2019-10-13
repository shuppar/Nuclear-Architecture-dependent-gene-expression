%% Removing nonunirom illuination

function IllCorrImage = hsvRet(I)

I = im2rgb(I);
hsv = rgb2hsv(I);
h = hsv(:, :, 1); % Hue image.
s = hsv(:, :, 2); % Saturation image.
v = hsv(:, :, 3); % Value (intensity) image

v = retinex_mccann99(v, 20);

hsv(:,:,3) = v;
IllCorrImage = hsv2rgb(hsv);
end