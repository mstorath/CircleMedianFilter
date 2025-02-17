function CMF_demoUnquantized

% load InSAR image, image source: https://earth.esa.int/workshops/ers97/papers/thiel/index-2.html
img = double(imread('../data/CMF_imgInSAR.png'));
% scale from [0, 255] to [-pi, pi]
img2pi = (img/255 - 0.5) * 2 * pi; 

% set filter size
R = 3;

% circle-median filter
tic
circleMedian = CMF_medfiltCirc2DMex(img2pi, R, R);
toc

% show results as hue value in hsv space
phaseToRGB = @(x) hsv2rgb(x/(2*pi) + 0.5, ones(size(x)), ones(size(x)));

subplot(1,2,1)
imshow(phaseToRGB(img2pi))
title('InSAR Data')

subplot(1,2,2)
imshow(phaseToRGB(circleMedian))
title('Arc distance median filter (3x3)')

end