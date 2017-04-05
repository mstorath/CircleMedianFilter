function phiAP = CMF_getAntipode( phi )
%CMF_getAntipode Computes the antipodal point of phi (-pi, pi]
phiAP = (phi - pi) .* (phi > 0) + (phi + pi) .* (phi <= 0);

end

