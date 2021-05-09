clear;

c = 343;
f_re = 2000;%частота для которой строиться МР
lamda = c/f_re;%длина волны 
N_one = 10;%кол-во датчиков
f_up = 4000;    %верхняя частота 
f_down = 0;  %нижняя частота
f_den = 0.5;    %шаг частот

d = lamda/2;%расстояние между решетками

thetaad = 0;       %пропускаемый угол
thetaan = 30;       %заглушаемый угол 

ula = phased.ULA(N_one,d);  %создаю равномерный линейный массив
ula.Element.BackBaffled = true;
w = [];     %массив весовых коэфф.
freq_range = f_down:f_den:f_up; %диапазон частот

%цикл для создания весовых коэффициентов для различных частот
for freq_range_while = f_down:f_den:f_up;    %прохожусь по каждой частоте
    lamda_sig = c/(freq_range_while);    %длина волны для определенной частоты
    
    wn = steervec(getElementPosition(ula)/lamda_sig,thetaan);   %вектор с пропускаемым углом

    wd = steervec(getElementPosition(ula)/lamda_sig,thetaad);   %вектор с подавляемым углом

    rn = wn'*wd/(wn'*wn);   %вычисление дроби из формулы

    wi = wd-wn*rn;          %вычисление самой формулы
    
    w = cat(1,w,wi');       %строю массив вес.коэффициентов
    
    
end


w_2zona = conj(w(8000:-1:2,:));
w(end,:) = real(w(end,:));

w_w = [w ; w_2zona];%весовые коэффициенты

%создаю ИХ
imp_har = ifft(w_w,length(w_w(:,1)),1);%ИХ всех весовых коэффициентов
for i = 1:10
    imp_har1(:,i) = fftshift(imp_har(:,i));%сдвигаю нулевой элемент в центр
end



%выделяю все ненулевые значения
n_down = 7950;
n_up = 8050;
imp_har_p = [];
for n_n = 1:10
    imp_har_p = [imp_har_p  imp_har1(n_down:1:n_up,n_n)];
    
end

for l = 1:10
    imp_har_p(:,l) = ifftshift(imp_har_p(:,l));%сдвигаю обратно нулевой элемент на место
end

fes_koef = fft(imp_har_p,length(imp_har_p(:,1)),1);%весовые коэффициенты

fes_koef_for_imp = fes_koef(1:51,:);
omega = -90:1:90;
ff_sig = 0:80:4000;

P = 0;
P_out = [];
%строю ДН
for j = 1:length(ff_sig)    %прохожусь по каждой частоте
    for i = 0:N_one-1   %считаю сумму дял ДН
    P_i = exp((-1i*i*2*pi*ff_sig(j)*d*(sind(omega)))/c)*(fes_koef_for_imp(j,i+1));%значение для каждого датчика по формуле (1.21)
    P = P + P_i;
    end
    BP = 20*log10((abs(P))/(max(abs(P))));%диагамма направлености по формуле (1.20)
    P_out = cat(1,P_out,BP);    %делаю массив для разных частот
    P = 0;
    BP = 0;
end 

f_up = 4000;
f_den = 80;
f_down = 0;
%%
[X,Y] = meshgrid(f_up:-f_den:f_down,omega);
C = X.*Y;
surf(X,Y,P_out.');
zlim([-60 0])
xlabel("Частота,Гц");
ylabel("Угол прихода, град");
zlabel("BP, дБ");
colormap('gray'); %палитра цветов
shading interp % убирает рёбра
grid on %включение сетки
colormap('gray'); %палитра цветов
caxis([-40, 0])
