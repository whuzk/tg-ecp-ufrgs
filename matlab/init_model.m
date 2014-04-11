[Fs,~,~,~,Data] = ecgutilities.interpret(EDB.e0103);
ECG = Data(1:5000,2);
N = length(ECG);
J = round(log2(Fs));

[lod,hid,lor,hir] = wfilters('bior3.5');
d0 = grpdelay(lod,1,1) + grpdelay(hid,1,1);
d1 = grpdelay(hid,1,1) + grpdelay(hir,1,1);
wdelay = round(2.^((J:-1:1)-1)*(d0+d1)-d0);
wdelay0 = wdelay(1);
wdelay(1) = [];