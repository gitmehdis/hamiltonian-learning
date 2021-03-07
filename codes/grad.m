function y = grad(beta,L,x)
[~,Obs] = DiagHamiltonian(beta,L,x);
y=-beta*Obs;
end
