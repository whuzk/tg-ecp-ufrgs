function J = pang_jpoints(R, RR, Fs)

% calcula a frequencia cardiaca e os intervalos
HR = RR*60./Fs;
A = HR < 100;
B = 100 <= HR & HR <= 110;
C = 100 <= HR & HR <= 11;
D = 120 < HR;

% calcula o ponto de medida de acordo com a frequencia cardiaca
L = zeros(size(R));
L(A) = 0.120;
L(B) = 0.112;
L(C) = 0.104;
L(D) = 0.100;

% calcula o ponto J
J = R + round(L.*Fs);