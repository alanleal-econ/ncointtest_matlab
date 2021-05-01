# ncointtest_matlab
ncointtest is a MATLAB implementation of the non-coointegration test, as developed by Souza et al (2018). The function takes five arguments:x and y, which are two time series, without missings, 
d, which is the order of integration of the two time series; b, which is the bandwidth in the test and it is a percentage of the number of observations, with n^b being the number of observations used in the test;
and r, which is the number of neighboring frequencies to be utilized in estimating the density spectral matrix.. 


Acknowledgements 

This command was translated from code belonging to Prof. Igor Viveiro Melo Souza, who graciouly allowed for its translation. MATLAB function for fractionary differentiation was based on
    freely distributed code in the fracdiff package in R. 

Disclaimer

    This program is provided without warranty of any kind. The author is not responsible for any cost derived by the usage of this program.


References

    Igor Viveiros Melo Souza, Valderio Anselmo Reisen, Glaura da Conceição Franco & Pascal Bondon (2018) The Estimation and Testing of the Cointegration Order Based on the Frequency
    Domain, Journal of Business & Economic Statistics, 36:4, 695-704, DOI: 10.1080/07350015.2016.1251442

    Leschinski, C., Voges, M. & Sibbertsen, P. (2020) A comparison of semiparametric tests for fractional cointegration. Stat Papers. https://doi.org/10.1007/s00362-020-01169-1


Author

    Alan Leal, University of São Paulo, Brazil
    prof@alanleal-econ.com
