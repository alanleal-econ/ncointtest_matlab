function [y_diff] = fracdiff_my(x_f,d)
% Fracdiff function, from fracdiff pacjkage in R. 

% More codes in: https://github.com/alanleal-econ and
% https://alanleal-econ.com

% The function has two arguments: x is a time series, both of a column vector and without missing
% values. d is the order of integration  to be considered in
% the calculation. The function returns a differenced series. 
%   Detailed explanation goes here

X_r=size(x_f,1);
x=x_f-mean(x_f);
PI=(1:X_r)';
PI(1)=-d;
for k=2:X_r
    PI(k)=PI(k-1)*(k-1-d)/k;
end
y_diff=x;

for i=2:X_r
   y_diff(i)=x(i)+sum(PI(1:(i-1)).*flip(x(1:(i-1))));
end
if (d~=0)
    y_diff=y_diff(2:X_r);
end
end

