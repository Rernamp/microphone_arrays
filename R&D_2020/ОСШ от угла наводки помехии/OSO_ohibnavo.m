clear;close all;

%заданию необходимые данные 

c = 343;
f_re = 2000;%частота для которой строиться МР
lamda = c/f_re;%длина волны 
N_one = 10;%кол-во датчиков
f_up = 4000;    %верхняя частота 
f_down = 0;  %нижняя частота
f_den = 0.5;    %шаг частот
freq_range = f_down:f_den:f_up; %диапазон частот
d = lamda/2;%расстояние между решетками
thetaad = 0;       %пропускаемый угол

osh_two = [];


    
[noise,fs] = audioread('nois.wav');%считываю шум

    t_duration = 2;%длительность сигнала
    t = 0:1/fs:t_duration-1/fs;
    NSampPerFrame = length(t);
    
hap = dsp.AudioPlayer('SampleRate',fs);
    %.................................Сигнал......................................
    dftFileReader = dsp.AudioFileReader('dft_voice_8kHz.wav',...
        'SamplesPerFrame',NSampPerFrame);
    sig = step(dftFileReader);
    

%созданию весовые коэф. и ИХ
ula = phased.ULA(N_one,d);  %создаю равномерный линейный массив
ula.Element.BackBaffled = true;
w = [];     %массив весовых коэфф.


thetaan = 30;
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

w_w = conj(w_w); %комплексно сопрягаю,в формуле так
%создаю ИХ
imp_har = ifft(w_w,length(w_w(:,1)),1);%ИХ всех весовых коэффициентов
for i = 1:10
    imp_har1(:,i) = fftshift(imp_har(:,i));    
end



%выделяю все ненулевые значения
n_down = 7950;
n_up = 8050;
imp_har_p = [];
for n_n = 1:10
    imp_har_p = [imp_har_p  imp_har1(n_down:1:n_up,n_n)];
    
end

for fi = 25:35

%принятие шума/создана что бы можно скрыть весь большой код
for p = 1
    N = 10;%кол-во решеток
    f_re = 2000;%частота
    %угол,под который пришел шум
    N_one = 10;%число отсчётов в одной последовательности

    %загружаю звуковые файлы
    
    

    

    nois = noise(1:length(sig));%укорачиваю кол-во отсчётов шума,что бы совпадали

    %определяю задержку
    tau = (d*sind(fi))/(c);

    dt = 1/fs;
    n_tau = (tau/dt); %перевожу в отсчёты



    spec = fft(nois);%БПФ шума

    %создаю массивы,которые дальше буду испольховать

    y_s_1zona = [];% первая зона Найквиста
    y_s_2zona = [];%вторая зона Найквиста
    nois_onMR = [];%сигнал на каждой решетке
    N_onezona =fix(length(sig)/(2*N_one))-1;%кол-во целых диапазонов в первой зоне Найквиста
    N_end1 = (fix(length(sig)/(2*N_one))*N_one)+2;%начало последнего диапазон
    N_end2 = (length(sig)/2)+1;%конец последнего диапазона
    N_centr = (N_end1 + N_end2)/2; % центральный отсчёт последнего диапазона

    %для каждой решетки смешаю на своё число отсчётов

    %приём шума
    for m = 1:N %цикл для каждой решетки
    
        for i = 0:N_onezona% цикл,который производит задержку шума
        
            n_tau_i = n_tau*(m-1);%кол-во задержаных отсчётов для каждой рещетки
            f_tau_i = ((2+2+N_one-1)/2 + i*N_one);%отсчёт,соответствующий центральной частоте для каждого диапазона
            k1 = spec(2+i*N_one:1+N_one+i*N_one) .* exp(-1i*2*pi*n_tau_i*f_tau_i/length(spec));%задержка/тут ---
            y_s_1zona = cat(1,y_s_1zona,k1);%присоединяю к пустому массиву

        end  

        y_s_1zona = [y_s_1zona (spec(N_end1:N_end2).*exp(-1i*2*pi*n_tau_i*N_centr/length(spec)))];%добавляю отсчёты,не вошедшие в основной массив

        y_s_1zona(end) = real(y_s_1zona(end));%так как отсчётов четное кол-во,то делаю реальным центральный отсчёт 

        y_s_2zona(1:N_end2 - 2) = conj(y_s_1zona((length(sig)/2)-1:-1:(length(sig)/2)-N_end2 + 2));%комплексно сопрягаю и отражаю отсчёты для второй зоны Найквиста

        y_f = [0 y_s_1zona.' y_s_2zona];%соединяю в отсчёты

        y_delayforMR = ifft(y_f); 

        y_delayforMR = y_delayforMR;%учитываю сигнал//что бы произвести дальнейшие действоя с звуковым сигналом,надо поменять y_delayforMR на sig.'       

        nois_onMR = cat(1,nois_onMR,y_delayforMR);%присоединяю полученный сигнал и так для каждой решетки



        %обнуляю массивы
        y_s_1zona = [];
        y_s_2zona = [];
        y_f = [];
        y_delayforMR = [];

    end

end

%убераю постоянку
nois_onMR = fft(nois_onMR,length(nois_onMR(1,:)),2);
nois_onMR(:,1) = 0;
nois_onMR = ifft(nois_onMR,length(nois_onMR(1,:)),2);

sig = fft(sig);
sig(1) = 0;
sig = ifft(sig);


%ОСШ без фильтрации



Power_one_exp_sig = mean(sig.^2);%мощность сигнала

nois_onMR_sre = mean(nois_onMR);%усредняю шум
Power_one_exp_nois = mean(nois_onMR_sre.^2);%мощность шума


OSH_one_exp_notfiltr = Power_one_exp_sig/Power_one_exp_nois;







%фильтрую сигнал
out_dec_sig = [];%сигнал на выходе
out_dec_nois = [];%шум на выходе

for m = 1:10
    out_d = filter(imp_har_p(:,m),1,nois_onMR(m,:));
    out_dec_nois = [out_dec_nois; out_d];
    out_d = filter(imp_har_p(:,m),1,sig.');
    out_dec_sig = [out_dec_sig; out_d];
    out_d = [];
    
end

%убераю постоянку
out_dec_sig = fft(out_dec_sig,length(out_dec_sig(1,:)),2);
out_dec_sig(:,1) = 0;
out_dec_sig = ifft(out_dec_sig,length(out_dec_sig(1,:)),2);


out_dec_nois = fft(out_dec_nois,length(out_dec_nois(1,:)),2);
out_dec_nois(:,1) = 0;
out_dec_nois = ifft(out_dec_nois,length(out_dec_nois(1,:)),2);

%усредняю сигнал и шум

out_dec_mean_sig = mean(out_dec_sig);

out_dec_mean_nois = mean(out_dec_nois);

%нахожу мощности после фильтрации 

Power_two_exp_sig = mean(out_dec_mean_sig.^2);

Power_two_exp_nois = mean(out_dec_mean_nois.^2);



OSH_two_exp_filtr = Power_two_exp_sig/Power_two_exp_nois;


osh_two = [osh_two OSH_two_exp_filtr];
end

osh_two = db(osh_two-OSH_one_exp_notfiltr);
angle = -5:5;
%%

plot(angle,osh_two);
grid on;
ylabel('Выигрыш в ОСШ, дБ');
xlabel("Ошибка наведения, град");
