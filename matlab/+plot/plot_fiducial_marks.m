function plot_fiducial_marks(data,F)

figure, plot(data);
hold on; grid on;

plot(F(1,:),data(F(1,:)),'sk');
plot(F(2,:),data(F(2,:)),'sr');
plot(F(3,:),data(F(3,:)),'ok');
plot(F(4,:),data(F(4,:)),'or');
plot(F(5,:),data(F(5,:)),'ok');
plot(F(6,:),data(F(6,:)),'^r');
plot(F(7,:),data(F(7,:)),'^k');

title('Fiducial marks');