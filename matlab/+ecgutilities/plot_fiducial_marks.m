function plot_fiducial_marks(data,R,F)

figure, plot(data);
hold on; grid on;

plot(R,data(R),'ok');
plot(F.P(:,1),data(F.P(:,1)),'sk');
plot(F.P(:,2),data(F.P(:,2)),'sk');
plot(F.R(:,1),data(F.R(:,1)),'ok');
plot(F.R(:,2),data(F.R(:,2)),'ok');
plot(F.T(:,1),data(F.T(:,1)),'^k');
plot(F.T(:,2),data(F.T(:,2)),'^k');
plot(F.IJ(:,1),data(F.IJ(:,1)),'xk');
plot(F.IJ(:,2),data(F.IJ(:,2)),'xk');

title('Fiducial marks');