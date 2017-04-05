function u = CMF_medfiltCircSep2D( y, R, T)
%CMF_medfiltCircSep2D Separable median filter for circle valued data based on
%projection of the componentwise median in R^2 

uCos = medfilt2(cos(y), [R, T], 'symmetric' );
uSin = medfilt2(sin(y), [R, T], 'symmetric' );
u = atan2(uSin, uCos);

% set data when atan2 undefined (should very rarely happen)
nDefIdx = find(max(abs(uCos), abs(uSin)) == 0);
u(nDefIdx) = y(nDefIdx);

end

