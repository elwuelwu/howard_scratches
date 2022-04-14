function Hn = hadamardGenNotNat(N)
m = log2(N);
if(N<2 || ~isequal(N, 2^nextpow2(N)))
    error('wrong input N!! N = 2^m for m=>1');
end
% a = 1:N;
% v = a;
Hn = 0;
for i=1:N
    aux1 = indexToStateMapper(i,m)';
    a    = de2bi(i-1,m,'left-msb')';
    for j=1:N
        v = de2bi(j-1,m,'left-msb');
        aux2 = indexToStateMapper(j,m);
%         Hn = Hn+((-1)^(mod(v*a,2)))*(mod(aux2*aux1,2));
        Hn = Hn+((-1)^(mod(v*a+v*v',2)))*(mod(aux2*aux1,2));
    end
end
Hn = Hn/sqrt(N);