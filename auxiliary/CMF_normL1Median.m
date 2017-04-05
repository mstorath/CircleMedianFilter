function med = CMF_normL1Median( y, varargin )
%CMF_normL1Median A wrapper for calculating the normalized L_1 median
% (which is based on the geometric median in R^2)

yEmbed = [cos(y(:)'); sin(y(:)')];
medEmbed = CMF_geometricMedianRN(yEmbed, varargin{:});
med = angle(medEmbed(1) + 1i*medEmbed(2));

end

