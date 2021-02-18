function ncointtest(x,y,d,band,r)
% Non-cointegration test as developed by Souza et al (2018)

% More codes in: https://github.com/alanleal-econ and
% https://alanleal-econ.com

% The function has five arguments: x and y are respectively the first and
% the second time series, both of a column vector and without missing
% values. d is the order of integration of both series to be considered in
% the test. b is the power of number of observations, a number between 0
% and 1, such as N^b is the number of observations used in the test. r is
% the number of neighbouring frequencies to be considered in the spectral
% density matrix. 

% Testing the input values:
% Testing for d
if (d<=0|d>1)
   disp("d must assume values between 0 and 1!");
   return;
end
% Testing for band
if (band<0|band>1)
   disp("Bandwidth must assume values between 0 and 1!");
   return;
end
% Testing for r:
if (r<=1)
   disp("r must be larger than 1!");
   return;
end

if ((sum(ismissing([x])+ismissing([y])))>=1)
   disp("Sample contains missing values!");
   return;
end


% 2) Differencing both series:
x=fracdiff_my(x,d);
y=fracdiff_my(y,d);

% 1) Implementing the test
x_r=size(x,1);
di11=(((1/(2*pi*x_r))*((fft(x).*conj(fft(x))))));di11=di11(2:end);
di12=(((1/(2*pi*x_r))*((fft(x).*conj(fft(y))))));di12=di12(2:end);
di21=(((1/(2*pi*x_r))*((fft(y).*conj(fft(x))))));di21=di21(2:end);
di22=(((1/(2*pi*x_r))*((fft(y).*conj(fft(y))))));di22=di22(2:end);
m=band;
band_temp=(floor(x_r^band))/r;
if (mod(band_temp,0.5)==0)
	band=floor(band_temp)*r;
else
	band=round(band_temp)*r;
end
di11_r=size(di11,1);
A=zeros(band,di11_r);
A_r_r=size(A,1)/r;

for i=1:A_r_r
   A(i,(r*i-r+1:r*i-r+1+r-1))=ones(1,r);
end

L_det=real(((A*di11).*(A*di22))-((A*di12).*(A*di21)));
L_det=log(L_det(1:A_r_r));
freq=(log(2.-2*cos(2*pi*((2:r:band)')/x_r)));
freq_r=size(freq,1);vec_1=ones(freq_r,1);
x_reg=[vec_1 freq];
y_reg=L_det;
b=((pinv(x_reg'*x_reg)*x_reg'*y_reg));b=b(2,1);
sd_as=sqrt((psi(1,3)+psi(1,2))./sum((freq-mean(freq)).^2));
t_ldr=b/sd_as;
p_value=1-normcdf(abs(t_ldr));

% 4) Printing the results of the test
disp("Non-cointegration test in the frequency domain");
disp("------------------------------------------------------------------------");
disp("H0: The bivariate time series vector is non-cointegrated");
disp("Note: Asymptotic values are based on Z ~ N(0,1)");
disp("------------------------------------------------------------------------");
fprintf("Set up: d = %g, Bandwidth = %g, r = %g, and n = %g \n",d,m,r,x_r);
disp("------------------------------------------------------------------------");
disp("        Est.b          ass.d        t-like stat      p-value")
disp("------------------------------------------------------------------------");
fprintf("%13.4g %15.4g %15.4g %15.4g \n",b,sd_as,t_ldr,p_value);
disp("------------------------------------------------------------------------");
end

