function Result = sogari_filter(Signal,Fs)

% linear filtering
Temp1 = smooth_and_diff(Signal,Fs);
%figure, plot(Temp1);

% non-linear transformation
Temp2 = mobd(Temp1,4);
%figure, plot(Temp2);

% shaping
Result = shape(Temp2,Fs);


function Result = smooth_and_diff(Signal,Fs)
Ws = round(0.01*Fs);
h1 = ones(1,Ws)/Ws;
h2 = wconv1(h1, [1 -1]);
Result = wconv1(Signal, h2, 'same');

function Result = mobd(Signal,M)
N = length(Signal);
Result = zeros(N,1);
for i = M:N
    w = Signal(i-M+1:i);
    p = abs(prod(w))*10^M;
    Result(i) = log(p+1)/log(M);
end

function Result = shape(Signal,Fs)
Ws = round(0.15*Fs);
h = ones(1,Ws);
Result = imdilate(Signal, h');
%Result = wconv1(Signal, h/Ws, 'same');
