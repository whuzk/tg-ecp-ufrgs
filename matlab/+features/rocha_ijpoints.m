function [I,J] = rocha_ijpoints(Beats, R, Fs)

N = 2^nextpow2(Fs/4);
Win = rectwin(N);
L = round(0.1*Fs);
M = size(Beats,2);
I = zeros(M,1);
J = zeros(M,1);
for i = 1:M
    % aplica a transformada de wigner-ville
    [WV,f] = math.wigner_ville(Beats(:,i), Win);
    S = sum(abs(WV(f < 0.2,:)));
    % extraçao do ponto isoeletrico
    [~,x] = min(S(R-L:R-1));
    I(i) = R - L + x - 1;
    % extraçao do ponto J
    [~,x] = min(S(R+1:R+L));
    J(i) = R + 1 + x - 1;
    % visualizacao
    %{
    figure(1), mesh(WV);
    figure(2), plot(WV(f < 0.2,:)');
    figure(3), plot(S); hold on;
    plot(I(i),S(I(i)),'ok');
    plot(J(i),S(J(i)),'ok');
    hold off;
    figure(4), plot(Beats(:,i)); hold on;
    plot(I(i),Beats(I(i),i),'ok');
    plot(J(i),Beats(J(i),i),'ok');
    hold off;
    pause;
    %}
end