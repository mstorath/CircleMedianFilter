function [fy, fx] = CMF_circMedFunc( y )
%CMF_circMedFunc computes the functional value of the defining functional
%of the circle-median for the data set y

% determine critical x-values, and sort them
y_unique = unique(y(:));
fx = sort([y_unique; CMF_getAntipode(y_unique)]);

% compute corresponding functional valued
fy = zeros(size(fx));
for i=1:numel(fx)
    fy(i) = sum(CMF_distCirc(fx(i), y));
end

end

