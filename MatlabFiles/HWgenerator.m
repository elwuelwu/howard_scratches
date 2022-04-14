function w = HWgenerator(a,b)
% this function produces D(a,b).
m = length(a);
if(m ~= length(b))
    error('Length of a and b mustb be equal!!');
end
x  = [0 1;1 0];
z  = [1 0;0 -1];
xz = x*z;
w = 1;
for i=1:m
    if(a(i) == 1 && b(i) == 1)
        w = kron(w,xz);        
    elseif(a(i) == 1 && b(i) == 0)
        w = kron(w,x);
    elseif(a(i) == 0 && b(i) == 1)
        w = kron(w,z);
    else
        w = kron(w,eye(2));
    end
end
