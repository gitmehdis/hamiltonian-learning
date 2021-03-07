function H=Hamiltonian(Jx,Jz,Jy,Jxx,Jzz,Jyy,L)
ld=2;
N = ld^L;
H = zeros(N,N);
Z=[1 0;0 -1];
X=[0 1;1 0];
Y=[0 -1i;1i 0];
XX=kron(X,X);
YY=kron(Y,Y);
ZZ=kron(Z,Z);


for i=1:L
    H=H + kron(eye(ld^(i-1)), kron( Jx(i)*X, eye(ld^(L-i))));
    H=H + kron(eye(ld^(i-1)), kron( Jz(i)*Z, eye(ld^(L-i))));
    H=H + kron(eye(ld^(i-1)), kron( Jy(i)*Y, eye(ld^(L-i))));
end

leftN = 1;
rightN = ld^(L-2);
for j=1:L-1
    h=(Jzz(j)*ZZ+Jxx(j)*XX+Jyy(j)*YY);
    H=H+kron(eye(leftN),kron(h,eye(rightN)));
    leftN = leftN*ld;
    rightN = rightN/ld;
end

H = H+Jzz(L).*kron(Z, kron(eye(ld^(L-2)), Z));
H = H+Jxx(L).* kron(X, kron(eye(ld^(L-2)), X));
H = H+Jyy(L).* kron(Y, kron(eye(ld^(L-2)), Y));

end
