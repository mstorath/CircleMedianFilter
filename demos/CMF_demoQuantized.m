function CMF_demoQuantized

% load data (data source: http://www.ndbc.noaa.gov/historical_data.shtml )
load CMF_windDir.mat

% convert from deg to rad in (-pi, pi]
myDegToRad = @(x)(x / 360 - 0.5) * 2 * pi;
windDir2Pi = myDegToRad(windDir);
% data is quantized to full degree
quant = myDegToRad(0:359);
% filter size around one day (data recorded every 10 minutes)
R = 24*6 + 1;
% perform filtering
tic
circleMedian = CMF_medfiltCirc2DQuantMex(windDir2Pi, R, 1, quant(:));
toc

% show result
subplot(2,1,1)
plot(windDir2Pi, '.')
title('Data')
axis tight

subplot(2,1,2)
plot(circleMedian, '.')
title('Arc distance median filter')
axis tight

end