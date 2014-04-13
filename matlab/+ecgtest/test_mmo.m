import ecgmmo.*
%close all;

% load ecg
[Fs,~,~,~,Data] = ecgutilities.interpret(EDB.e0103);
Signal = Data(1:1000,1);

% apply operators
s = 10;
B = zeros(1,2*s+1);
Y1 = mmo_erosion(Signal,B);
Y2 = mmo_dilation(Signal,B);
Y3 = mmo_opening(Signal,B,B);
Y4 = mmo_closing(Signal,B,B);
Y5 = mmo_derivative(Signal,B);
Y6 = mmo_derivative_flat(Signal,s);
Y7 = mmo_open_closing(Signal,B,B);
Y8 = mmo_close_opening(Signal,B,B);
Y9 = mmo_average(Signal,B,B);
Y10 = rterosion(Signal,B);
Y11 = rtdilation(Signal,B);
Y12 = rtopening(Signal,B,B);
Y13 = rtclosing(Signal,B,B);
Y14 = rtopenclosing(Signal,B,B);
Y15 = rtcloseopening(Signal,B,B);
Y16 = rtaverage_flat(Signal,s,s);

% plot figures
figure, plot(Signal);
figure, plot([Y1 Y2]);
figure, plot([Y3 Y4]);
figure, plot([Y5 Y6]);
figure, plot([Y7 Y8]);