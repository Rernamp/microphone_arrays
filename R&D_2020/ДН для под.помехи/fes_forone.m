clear;

f_re = 2000;%частота дл€ которой строитьс€ ћ–
c = 343;%скорость света
lamda = c/f_re;%длина волны 
N_one = 10;%кол-во датчиков
omega = -90:2:90;
f_up = 4000;    %верхн€€ частота 
f_down = 300;  %нижн€€ частота
f_den = 100;    %шаг частот

d = lamda/2;%рассто€ние между решетками

thetaad = 10;       %пропускаемый угол
thetaan = 40;       %заглушаемый угол 

ula = phased.ULA(N_one,d);  %создаю равномерный линейный массив
ula.Element.BackBaffled = true;
w = [];     %массив весовых коэфф.
ff_sig = f_down:f_den:f_up; %диапазон частот

%цикл дл€ создани€ весовых коэффициентов дл€ различных частот
for fff_sig = f_down:f_den:f_up;    %прохожусь по каждой частоте
    lamda_sig = c/(fff_sig);    %длина волны дл€ определенной частоты
    
    wn = steervec(getElementPosition(ula)/lamda_sig,thetaan);   %вектор с пропускаемым углом

    wd = steervec(getElementPosition(ula)/lamda_sig,thetaad);   %вектор с подавл€емым углом

    rn = wn'*wd/(wn'*wn);   %вычисление дроби из формулы

    wi = wd-wn*rn;          %вычисление самой формулы
    
    w = cat(1,w,wi');       %строю массив вес.коэффициентов
    
end



%вспомогательные переменные
P = 0;
P_out = [];

%вычисл€ю диаграмму направлености
for j = 1:length(ff_sig)    %прохожусь по каждой частоте
    for i = 0:N_one-1   %считаю сумму д€л ƒЌ
    P_i = exp((-1i*i*2*pi*ff_sig(j)*d*(sind(omega)))/c)*w(j,i+1);%значение дл€ каждого датчика по формуле (1.21)
    P = P + P_i;
    end
    BP = 20*log10((abs(P))/(max(abs(P))));%диагамма направлености по формуле (1.20)
    P_out = cat(1,P_out,BP);    %делаю массив дл€ разных частот
    P = 0;
    BP = 0;
end 


%%

figure(1);
hold on;
plot(90:-2:-90,P_out(1,:));
plot(90:-2:-90,P_out(2,:));
plot(90:-2:-90,P_out(3,:));
plot(90:-2:-90,P_out(4,:));
plot(90:-2:-90,P_out(5,:));
plot(90:-2:-90,P_out(6,:));
plot(90:-2:-90,P_out(7,:));
plot(90:-2:-90,P_out(8,:));
plot(90:-2:-90,P_out(9,:));
plot(90:-2:-90,P_out(10,:));
plot(90:-2:-90,P_out(11,:));
plot(90:-2:-90,P_out(12,:));
plot(90:-2:-90,P_out(13,:));
plot(90:-2:-90,P_out(14,:));
plot(90:-2:-90,P_out(15,:));
plot(90:-2:-90,P_out(16,:));
plot(90:-2:-90,P_out(17,:));
plot(90:-2:-90,P_out(18,:));
plot(90:-2:-90,P_out(19,:));
plot(90:-2:-90,P_out(20,:));
plot(90:-2:-90,P_out(21,:));
plot(90:-2:-90,P_out(22,:));
plot(90:-2:-90,P_out(23,:));
plot(90:-2:-90,P_out(24,:));
plot(90:-2:-90,P_out(25,:));
plot(90:-2:-90,P_out(26,:));
line([10,10], [-400,30]);
line([40,40], [-400,30]);
ylim([-300,30]);
grid on;
xlim([-90,90]);
ylim([-100, 10])
xlabel("”гол прихода, град");
ylabel("BP, дЅ");

figure(2);
[X,Y] = meshgrid(f_up:-f_den:f_down,90:-2:-90);
surf(X,Y,P_out.');
xlabel("„астота,√ц");
ylabel("”гол прихода, град");
zlabel("BP, дЅ");
shading interp; % убирает рЄбра
ylim([-90,90]);
zlim([-40,0])


grid on %включение сетки
colormap('gray'); %палитра цветов
caxis([-40, 0])

