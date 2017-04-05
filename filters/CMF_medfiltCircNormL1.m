function u = CMF_medfiltCircNormL1( y, R, T, maxIter, stopTol)
%CMF_medfiltCircNormL1 Computes the normalized L_1 median filter for circle valued data. 
% It is based on projection of the geometric median in R^2 (using Weiszfeld algorithm)

% default number iterations
if ~exist('maxIter', 'var')
    maxIter = 100;
end

% default stopping criterion
if ~exist('stopTol', 'var')
    stopTol = 1e-9;
end

% compute geometric median in R^2 (using fast C++ implementation)
yR2 = cat(3, cos(y), sin(y));
geoMedR2 = CMF_medfiltGeoRN2DMex(yR2, R, T, maxIter, stopTol);

% projection to circle
u = angle(geoMedR2(:,:,1) + 1i*geoMedR2(:,:,2));

% set data when atan2 undefined (should very rarely happen)
nDefIdx = find(max(abs(geoMedR2), [], 3) == 0);
u(nDefIdx) = y(nDefIdx);

end

