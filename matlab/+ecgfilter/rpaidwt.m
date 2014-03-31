function a = rpaidwt(a,D,N,H,G)

J = length(D);
A = cell(J+1,1);
A{J+1} = a;
D = [0; D];
H = H(:);
G = G(:);
A = ridwt(0,J-1,J,A,D,H,G);
for i = 1:N-2
    A = ridwt(i,1,J,A,D,H,G);
end
a = A{1}(:);

function A = ridwt(i,j,J,A,D,H,G)
if j >= J
    return;
elseif (i == 0)
    if (j == 0)
        A{1+0}(1+0) = A{1+1}(1+0)*H(1) + D{1+1}(1+0)*G(1);
        A{1+0}(1+1) = A{1+1}(1+0)*H(2) + D{1+1}(1+0)*G(2);
    else
        A{1+j}(1+0) = A{1+j+1}(1+0)*H(1) + D{1+j+1}(1+0)*G(1);
        A = ridwt(0,j-1,J,A,D,H,G);
    end
elseif mod(i,2) ~= 0
    k = (i+1)/2;
    if mod(k,2) == 0
        p = k/2;
        m = 1:min(p,ceil(length(H)/2));
        A{1+j}(1+k) = A{1+j+1}(p-m+1)*H(2*m-1) + D{1+j+1}(p-m+1)*G(2*m-1);
    else
        p = (k-1)/2;
        m = 1:min(p,floor(length(H)/2));
        A{1+j}(1+k) = A{1+j+1}(p-m+1)*H(2*m) + D{1+j+1}(p-m+1)*G(2*m);
    end
    if j == 1
        m = 1:min(k,ceil(length(H)/2));
        A{1+0}(1+i+1) = A{1+1}(k-m+1)*H(2*m-1) + D{1+1}(k-m+1)*G(2*m-1);
        m = 1:min(k,floor(length(H)/2));
        A{1+0}(1+i+2) = A{1+1}(k-m+1)*H(2*m) + D{1+1}(k-m+1)*G(2*m);
    end
else
    A = ridwt(i/2,j+1,J,A,D,H,G);
end