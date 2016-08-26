function coeff = barycoefs( offset, kappa, dist)
% BARYCOEFS - compute local Lagrange interpolation
%
% coeff = barycoefs(offset,kappa,index) returns coeff given the offset,
% target point kappa, and the source from which we are taking the function
% value.
% Hai 05/24/16.
%   
n = numel(offset);
whts = ones(n,1);

for i = 1:n
    for j = 1:n
        if i~=j
            whts(i) = whts(i)/(offset(i)-offset(j));
        end
    end
end

dd = 0*kappa;
for i = 1:n
    dd = dd + whts(i)./(kappa-offset(i));
end
index = dist + ceil(n/2);
coeff = (whts(index)./(kappa-offset(index))) ./dd;

end

