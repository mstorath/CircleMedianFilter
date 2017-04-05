function d = CMF_distCirc( phi, psi )
%CMF_distCirc Computes the circular distance of phi and psi 
d = abs(phi - psi);
d(d > pi) = 2 * pi - d(d > pi);

end

