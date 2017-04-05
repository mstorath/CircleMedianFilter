function med = CMF_separableMedian( y )
%CMF_separableMedian Computes the separable median;
%that is, copmuting the median of each component separately and
%and projecting the result to the unit circle

med1 = median(cos(y(:)));
med2 = median(sin(y(:)));
med = angle(med1 + 1i*med2);

end

