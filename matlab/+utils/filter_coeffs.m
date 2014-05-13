function Result = filter_coeffs

bs50 = zeros(20,5);
as50 = zeros(20,5);
ds50 = zeros(20,1);
bl50 = zeros(20,5);
al50 = zeros(20,5);
dl50 = zeros(20,1);
Fm = 50;
for i = 3:20
    Fs = i*Fm;
    [bs50(i,:),as50(i,:)] = cheby1(2,0.1,[Fm-1 Fm+1]*2/Fs,'stop');
    ds50(i) = mean(grpdelay(bs50(i,:),as50(i,:),[5 15]*2/Fs));
    [bl50(i,:),al50(i,:)] = butter(4,40*2/Fs);
    dl50(i) = mean(grpdelay(bl50(i,:),al50(i,:),[5 15]*2/Fs));
end

bs60 = zeros(16,5);
as60 = zeros(16,5);
ds60 = zeros(16,1);
bl60 = zeros(16,5);
al60 = zeros(16,5);
dl60 = zeros(16,1);
Fm = 60;
for i = 3:16
    Fs = i*Fm;
    [bs60(i,:),as60(i,:)] = cheby1(2,0.1,[Fm-1 Fm+1]*2/Fs,'stop');
    ds60(i) = mean(grpdelay(bs60(i,:),as60(i,:),[5 15]*2/Fs));
    [bl60(i,:),al60(i,:)] = butter(4,40*2/Fs);
    dl60(i) = mean(grpdelay(bl60(i,:),al60(i,:),[5 15]*2/Fs));
end

Result.bs50 = bs50;
Result.bs50 = as50;
Result.bs50 = ds50;
Result.bs50 = bl50;
Result.bs50 = al50;
Result.bs50 = dl50;

Result.bs50 = bs60;
Result.bs50 = as60;
Result.bs50 = ds60;
Result.bs50 = bl60;
Result.bs50 = al60;
Result.bs50 = dl60;