function x = CMF_wrapAngle( y )
%CMF_wrapAngle Wraps an angle to the interval [-pi, pi]
x = angle(exp(1i * y));

end

