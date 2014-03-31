function [A,D] = rtdwt(X,J,H,G,Lo_R,Hi_R)

N = length(X);
L = length(H);
A = cell(J,1);
D = cell(J,1);
S = zeros(J+1,L);
len = N;
for j = 1:J
    len = floor(len/2);
    A{j} = zeros(1,len);
    D{j} = zeros(1,len);
end

%{
Signal1 = zeros(1,N);
Signal2 = zeros(1,N);
figure(1);
[lh1,lh2,dwth] = setup_plot(Signal1,Signal2,D,J,N);
%}

X = [X; zeros(200,1)];
k = 1;
for i = 1:length(X)
    S(1,:) = push(S(1,:), X(i));
    j = levelof(k);
    if j <= J
        [S,A,D] = update_dwt(j,S,A,D,H,G);
        k = k + 1;
    else
        %{
        DD = D;%ecgfilter.denoise(D,J,N);
        Signal2 = ecgfilter.dpaidwt(A{end},DD,N,Lo_R,Hi_R);
        set(lh2, 'ydata', Signal2);
        ecgutilities.plot_dwt(DD,[],J,N,dwth);
        pause;
        %}
        k = 1;
    end
    %{
    Signal1 = push(Signal1, X(i));
    set(lh1, 'ydata', Signal1);
    drawnow;
    %}
end


function [S,A,D] = update_dwt(j,S,A,D,H,G)
a = S(j,:)*H(end:-1:1)';
d = S(j,:)*G(end:-1:1)';
if j < size(S,1)
    S(j+1,:) = push(S(j+1,:), a);
end
A{j} = push(A{j}, a);
D{j} = push(D{j}, d);

function j = levelof(k)
j = 1;
while mod(k,2) == 0
    j = j + 1;
    k = k / 2;
end

function Y = push(X,v)
Y = [X(length(v)+1:end) v];

function [lh1,lh2,dwth] = setup_plot(Signal1,Signal2,D,J,N)
subplot(3,1,1); grid on;
title('Original signal');
lh1 = line((1:N)',Signal1);
xlim([1 N]);
subplot(3,1,2); grid on;
title('Reconstructed signal');
lh2 = line((1:N)',Signal2);
xlim([1 N]);
subplot(3,1,3);
dwth = ecgutilities.plot_dwt(D,[],J,N,[]);