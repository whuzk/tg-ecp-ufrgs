function Result = extract_features(MethodName, Data)

switch MethodName
    case 'Rocha'
        Result = rocha_features(Data.Beats, Data.Fs, Data.RR);
    case 'Mohebbi'
        Result = mohebbi_features(Data.Beats, Data.Fs, Data.Template);
    case 'Gopalak'
        Result = gopalak_features(Data.Beats);
    otherwise
        error('invalid method name');
end


function Result = rocha_features(Beats, Fs, RR)
I = [];
J = [];
[B1,B2] = ecgfilter.st_deviation(Beats, Fs, RR);
[C1,C2] = ecgfilter.hermite_coeff(Beats, I, J, min(RR,Fs));
Result = [B1 B2 C1 C2];

function Result = mohebbi_features(Beats, Fs, Template)
STlen = fix(0.16*Fs);
Jtemp = ecgfilter.jay_point(Template, Fs);
STtemp = Template(Jtemp:Jtemp+STlen-1);
modelF = (STtemp(1:2:end) + STtemp(2:2:end))/2;

n = length(modelF);
m = size(Beats,2);
Result = zeros(m,n);
for i = 1:m
   J = ecgfilter.jay_point(Beats(:,i), Fs);
   ST = Beats(J:J+STlen-1,i);
   F = (ST(1:2:end) + ST(2:2:end))/2;
   Result(i,:) = modelF - F;
end

function Result = gopalak_features(Beats)
n = size(Beats,1);
H = ecgmath.hermite(n,50,1.5);
Result = Beats'*H;
