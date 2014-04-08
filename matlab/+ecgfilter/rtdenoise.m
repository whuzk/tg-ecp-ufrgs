function Result = rtdenoise(X,J,M,W)

Lo_R = sqrt(2)*W;
Hi_R = qmf(Lo_R);
Hi_D = wrev(Hi_R);
Lo_D = wrev(Lo_R);

N = length(X);
Result = zeros(N,1);
X = [zeros(1,M/2-1) X(:)' zeros(1,M/2)];
for i = 1:N
    [A,D] = ecgfilter.dpadwt1(X(i:i+M-1),J,Lo_D,Hi_D);
    D = ecgfilter.denoise(D,J,M);
    R = ecgfilter.dpaidwt1(A{end},D,M,Lo_R,Hi_R);
    Result(i) = R(end/2);
end