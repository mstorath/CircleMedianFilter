function med = CMF_arcdistMedian( y )
%CMF_arcdistMedian 
% Brute force method for computing an arc distance median of y

n = length(y);
dCirc = zeros(n,1);
for i = 1:numel(y)
    dCirc(i) = sum(CMF_distCirc(y(i),y));
end
[~, minIdx] = min(dCirc);
med = y(minIdx);

end

