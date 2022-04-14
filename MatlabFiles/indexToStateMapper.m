function v = indexToStateMapper(index,m)
if(index > 2^m)
    error('index must be less than 2^m');
end
v0   = [1;0];
v1   = [0;1];
auxv = de2bi(index-1,m, 'left-msb');
v    = 1;
for i=1:m
    if(auxv(i) == 0)
        v = kron(v,v0);
    else
        v = kron(v,v1);
    end
end