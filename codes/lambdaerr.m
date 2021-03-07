function Obserr = lambdaerr(beta,L,coeff,errparam)
[~,Obs] = DiagHamiltonian(beta,L,coeff);
noise=randn(6*L,1).*errparam;
Obserr=noise+Obs;
