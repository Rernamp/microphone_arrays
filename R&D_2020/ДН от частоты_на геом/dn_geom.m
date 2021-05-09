clear;

f_re = 3000;%частота дл€ которой строитьс€ ћ–
c = 343;%скорость света
lamda = c/f_re;%длина волны 
N_one = 10;%кол-во датчиков
omega = -90:2:90;
f_up = 4000;    %верхн€€ частота 
f_down = 100;  %нижн€€ частота
f_den = 100;    %шаг частот

d = lamda/2;%рассто€ние между решетками


w = [];     %массив весовых коэфф.
ff_sig = f_down:f_den:f_up; %диапазон частот
ww = [];
ksi = 0;
%цикл дл€ создани€ весовых коэффициентов дл€ различных частот
    %прохожусь по каждой частоте
    
       %длина волны дл€ определенной частоты
F = zeros(1,length(omega));
Ff = [];
for fff_sig = f_down:f_den:f_up
    lamda_i = c/fff_sig;
    for i = 1:N_one
        F = F + ((-exp((-1i*2*pi*d*(i-1)*sind(omega-ksi))/(lamda_i)))/((N_one)^0.5));
    end
    Ff = [Ff ; F];
    F = F*0;
end



set(0,'DefaultAxesFontSize',14,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',14,'DefaultTextFontName','Times New Roman'); 
kF = abs(Ff);
for fff_sig = (f_down:f_den:f_up)/f_down
    kF(fff_sig,:) = kF(fff_sig,:)/max(kF(fff_sig,:));
end
    
   
%%

[X,Y] = meshgrid(f_up:-f_den:f_down,omega);
C = X.*Y;
surf(X,Y,db(kF).');
zlim([-40 0]);
xlabel("„астота,√ц");
ylabel("”гол прихода, град");
colormap('gray'); %палитра цветов
caxis([-40, 0])
zlabel("BP, дЅ");
colormap('gray'); %палитра цветов
shading interp % убирает рЄбра
grid on %включение сетки
