function nric = CMF_nric( data, groundTruth, estimate )
%CMF_nric Computes the circular noise reduction index (NRI_c)
%
% See Nikolaidis, Pitas, "Nonlinear Processing and Analysis of Angular
% Signals", 1998

num = sum( CMF_distCirc(groundTruth(:), data(:)) );
den = sum( CMF_distCirc(groundTruth(:), estimate(:))  );
nric = 10 * log10 ( num/den );

end

