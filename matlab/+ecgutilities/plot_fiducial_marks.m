function plot_fiducial_marks(data,F)

figure, plot(data);
hold on; grid on;
plot(F.P(1),data(F.P(1)),'sk');
plot(F.P(2),data(F.P(2)),'sk');
plot(F.P(3),data(F.P(3)),'sk');
plot(F.R(1),data(F.R(1)),'ok');
plot(F.R(2),data(F.R(2)),'ok');
plot(F.R(3),data(F.R(3)),'ok');
plot(F.T(1),data(F.T(1)),'^k');
plot(F.T(2),data(F.T(2)),'^k');
plot(F.T(3),data(F.T(3)),'^k');