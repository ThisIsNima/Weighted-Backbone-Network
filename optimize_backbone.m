function F = optimize_backbone(std, std_vals) 

%This is the optimizationfunction that is solved within the system of N equations 

amat = std*std';

MLmat=amat-std_vals;

%MLmat = (tau.*amat)./(ones(Nd,Nd) - amat);

MLmat(logical(eye(size(MLmat)))) = nan;% diagonal be nan.

F = sum(MLmat,2,'omitnan');
