function u = Tridiag(a,b,c,r)
%
% Resolution of a tridiagonal sustem with matrix
%
%  | b_1 c_1 0   ...                  |
%  | a_2 b_2 c_2 ...                  |
%  |             ...                  |
%  |             ... a_N-1 b_N-1 c_N-1|
%  |             ...    0  a_N   b_N  |
%
% Courtesy of Mattia Zanella (University of Pavia).

N=length(b);
bet=b(1);
u(1)=r(1)/b(1);
for j=2:N
   g(j)=c(j-1)/bet;
   bet=b(j)-a(j)*g(j);
   u(j)=(r(j)-a(j)*u(j-1))/bet;
end
for j=N-1:-1:1
   u(j)=u(j)-g(j+1)*u(j+1);
end