function pnt = CMF_angleToPointR2( ang )
%CMF_ANGLETOPOINTR2 Converts an angle to its representation as unit vector 
%in R2, saved as a numel(ang) x 2 array
pnt = [cos(ang(:)'); sin(ang(:)')];

end

