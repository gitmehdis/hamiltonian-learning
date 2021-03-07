function [logZ,Obs]=DiagHamiltonian(beta,L,coeff)

ld=2;
% Recall that the boundary conditions are set in the coefficients below

Jxx=coeff(1:L,1);
Jyy=coeff(L+1:2*L,1);
Jzz=coeff(2*L+1:3*L,1);
Jx=coeff(3*L+1:4*L,1);
Jy=coeff(4*L+1:5*L,1);
Jz=coeff(5*L+1:6*L,1);
 
N = ld^L;
H = zeros(N,N);
Z=[1 0;0 -1];
X=[0 1;1 0];
Y=[0 -1i;1i 0];
XX=kron(X,X);
YY=kron(Y,Y);
ZZ=kron(Z,Z);


H=Hamiltonian(Jx,Jz,Jy,Jxx,Jzz,Jyy,L);
[eigvec,en] = eig(H);
en=diag(en);

invV = 1.0/L;
rho=zeros(ld^L);
partfun=0;
for k=1:size(en)
    coeff = exp(-beta*en(k));
    partfun = partfun+coeff;
    rho = rho + coeff.*eigvec(:,k)*(eigvec(:,k)');
end
invZ = 1.0/partfun;
logZ = log(partfun);



ex=zeros(L,1);
ey=zeros(L,1);
ez=zeros(L,1);
lN = 1;
rN = ld^(L-1);
for j=1:L
    ObsZ = kron(eye(lN), kron(Z, eye(rN)));
    ObsX = kron(eye(lN), kron(X, eye(rN)));
    ObsY = kron(eye(lN), kron(Y, eye(rN)));
    ez(j)=trace(ObsZ*rho)*invZ;
    ex(j)=trace(ObsX*rho)*invZ;
    ey(j)=trace(ObsY*rho)*invZ;
    lN = lN*ld;
    rN = rN/ld;
end

lN = 1;
rN = ld^(L-2);
ezz = zeros(L,1);
exx = zeros(L,1);
eyy = zeros(L,1);
for j=1:L-1
    ObsZ = kron(eye(lN), kron(kron(Z,Z), eye(rN)));
    ObsX = kron(eye(lN), kron(kron(X,X), eye(rN)));
    ObsY = kron(eye(lN), kron(kron(Y,Y), eye(rN)));
    ezz(j)=trace(ObsZ*rho)*invZ;
    exx(j)=trace(ObsX*rho)*invZ;
    eyy(j)=trace(ObsY*rho)*invZ;
    lN = lN*ld;
    rN = rN/ld;
end
ObsZ = kron(Z, kron(eye(ld^(L-2)),Z));
ObsX = kron(X, kron(eye(ld^(L-2)),X));
ObsY = kron(Y, kron(eye(ld^(L-2)),Y));

ezz(L) = trace(ObsZ*rho)*invZ;
exx(L) = trace(ObsX*rho)*invZ;
eyy(L) = trace(ObsY*rho)*invZ;


Obs=[exx;eyy;ezz;ex;ey;ez];
% [logZ exx, eyy, ezz]
end
