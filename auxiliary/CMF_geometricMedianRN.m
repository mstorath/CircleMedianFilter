function med = CMF_geometricMedianRN( y, nIter )
%CMF_geometricMedianRN A simple implementation of the Weiszfeld algorithm
%for the geometric median in R^M
% y in Dimension M x Samples N

% default number of iterations
if ~exist('nIter', 'var')
   nIter = 100;
end

% init
[M, N] = size(y);
med = zeros(M, 1);

% Weiszfeld iteration
for i = 1:nIter
   num = zeros(M, 1);
   den = 0;
   for n = 1:N
       dist = norm(y(:,n) - med, 2);
       num = num + y(:,n) / dist;
       den = den + 1 / dist;
   end
   if ~isfinite(den)
       break;
   end
   med = num/den;
end
    

end

