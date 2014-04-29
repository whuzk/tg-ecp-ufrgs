function plot_fiducial_marks(data,F)

figure, plot(data);
hold on; grid on;

for i = 1:length(F)
    plot(F(i).P(1),data(F(i).P(1)),'sk');
    plot(F(i).P(2),data(F(i).P(2)),'sk');
    plot(F(i).P(3),data(F(i).P(3)),'sk');
    plot(F(i).R(1),data(F(i).R(1)),'ok');
    plot(F(i).R(2),data(F(i).R(2)),'ok');
    plot(F(i).R(3),data(F(i).R(3)),'ok');
    plot(F(i).T(1),data(F(i).T(1)),'^k');
    plot(F(i).T(2),data(F(i).T(2)),'^k');
    plot(F(i).T(3),data(F(i).T(3)),'^k');
end