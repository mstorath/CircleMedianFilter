function medSet = CMF_bisectMedian( y, testSetAdd )
%CMF_bisectMedian Computes a set of bisecting medians of y
% y: data
% testSetAdd: additional set of candidates to be tested for being a bisecting median, default: []

medSet = [];
n = length(y);
n2 = (n/2);

if ~exist('testSetAdd', 'var')
   testSetAdd = []; 
end
% test set consists of points in y, its antipodes, and an additional test
% set of values to be tested
testSet = [y(:); CMF_getAntipode(y(:)); testSetAdd(:)];

% checks which points in the test set are a bisecting median
m = numel(testSet);
for i=1:m
    candPoint = testSet(i);
    a = angle( exp(1i * (candPoint - y)));
    cclw = sum((0 <= a)&(a < pi)) + sum(a == pi) + sum(a == -pi);
    clw = sum((-pi < a)&(a <= 0)) + sum(a == pi) + sum(a == -pi);
    if (cclw >= n2) && (clw >= n2)
       
       numThis = sum(CMF_distCirc(candPoint, y) <= pi/2);
       numAP = sum(CMF_distCirc(candPoint, y) >= pi/2);
       if numThis >= numAP
           medSet = [medSet; candPoint];
       end
    end
end


end

