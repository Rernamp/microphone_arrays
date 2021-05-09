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

wn = steervec(getElementPosition(ula)/lamda_sig,thetaan)
w_2zona = conj(w(8000:-1:2,:));
w(end,:) = real(w(end,:));

w_w = [w ; w_2zona];%весовые коэффициенты

%создаю ИХ
imp_har = ifft(w_w,length(w_w(:,1)),1);%ИХ всех весовых коэффициентов
for i = 1:10
    imp_har1(:,i) = fftshift(imp_har(:,i));
end
imp_har_dec = [];
imp_har_manual = [];
imp_har_dontdec = [];
N_imp_har = 32;
range_one = (length(w_w(:,1)))/N_imp_har;
%создаю ИХ для каждого элемента МР
for i = 1:10
    imp_h = decimate(imp_har(:,i),range_one);%ИХ через готовую функцию
    imp_har_dec = [imp_har_dec imp_h];
    imp_k = imp_har(range_one:range_one:16000,i);%ИХ при помощи ручной выборки
    imp_har_manual = [imp_har_manual imp_k];
    imp_har_l = imp_har(1:N_imp_har,i);%ИХ выбираю первый 32 элемента
    imp_har_dontdec = [imp_har_dontdec imp_har_l];
    imp_h = [];
    imp_k = [];
    imp_har_l = [];
end

%принятие сигнала/создана что бы можно скрыть весь большой код
for p = 1
    N = 10;%кол-во решеток
    f_re = 2000;%частота
    fi = 30;%угол,под который пришел шум
    N_one = 10;%число отсчётов в одной последовательности

    %загружаю звуковые файлы
    
    [noise,fs] = audioread('nois.wav');%считываю шум

    t_duration = 2;%длительность сигнала
    t = 0:1/fs:t_duration-1/fs;
    NSampPerFrame = length(t);

    hap = dsp.AudioPlayer('SampleRate',fs);
    %.................................Сигнал......................................
    dftFileReader = dsp.AudioFileReader('dft_voice_8kHz.wav',...
        'SamplesPerFrame',NSampPerFrame);
    sig = step(dftFileReader);

    nois = noise(1:length(sig));%укорачиваю кол-во отсчётов шума,что бы совпадали

    %определяю задержку
    tau = (d*sind(fi))/(c);

    dt = 1/fs;
    n_tau = (tau/dt); %перевожу в отсчёты



    spec = fft(nois);%БПФ шума

    %создаю массивы,которые дальше буду испольховать

    y_s_1zona = [];% первая зона Найквиста
    y_s_2zona = [];%вторая зона Найквиста
    y_s_oneMR = [];%сигнал на каждой решетке
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
            k1 = spec(2+i*N_one:1+N_one+i*N_one) .* exp(-1i*2*pi*n_tau_i*f_tau_i/length(spec));%задержка
            y_s_1zona = cat(1,y_s_1zona,k1);%присоединяю к пустому массиву

        end  

        y_s_1zona = [y_s_1zona (spec(N_end1:N_end2).*exp(-1i*2*pi*n_tau_i*N_centr/length(spec)))];%добавляю отсчёты,не вошедшие в основной массив

        y_s_1zona(end) = real(y_s_1zona(end));%так как отсчётов четное кол-во,то делаю реальным центральный отсчёт 

        y_s_2zona(1:N_end2 - 2) = conj(y_s_1zona((length(sig)/2)-1:-1:(length(sig)/2)-N_end2 + 2));%комплексно сопрягаю и отражаю отсчёты для второй зоны Найквиста

        y_f = [spec(1) y_s_1zona.' y_s_2zona];%соединяю в отсчёты

        y_delayforMR = ifft(y_f); 

        y_delayforMR = y_delayforMR + sig.';%учитываю сигнал//что бы произвести дальнейшие действоя с звуковым сигналом,надо поменять y_delayforMR на sig.'       

        y_s_oneMR = cat(1,y_s_oneMR,y_delayforMR);%присоединяю полученный сигнал и так для каждой решетки



        %обнуляю массивы
        y_s_1zona = [];
        y_s_2zona = [];
        y_f = [];
        y_delayforMR = [];

    end

end

%фильтрую сигнал
out_dec = [];
out_dontdec = [];
out_manual = [];
for m = 1:10
    out_d = filter(imp_har_dec(:,m),1,y_s_oneMR(m,:));
    out_dec = [out_dec; out_d];
    out_dd = filter(imp_har_dontdec(:,m),1,y_s_oneMR(m,:));
    out_dontdec = [out_dontdec; out_dd];
    out_m = filter(imp_har_manual(:,m),1,y_s_oneMR(m,:));
    out_manual = [out_manual; out_m];
end
%нахожу среднее
out_dec_mean = mean(out_dec);
out_dontdec_mean = mean(out_dontdec);
out_manual_mean = mean(out_manual);
% %строю профильтрованые сигналы
% figure(1);
% title('Сигналы после фильтрации');
% hold on;
% plot(t,out_dec_mean);
% plot(t,out_dontdec_mean);
% plot(t,out_manual_mean);
% legend('децимация через функцию','децимация вручную','первые 32 отсчёта');
% hold off;
%%

tt = 1:80;
figure(2);

plot(tt,imp_har(tt,:));
grid on;
ylabel("Амплитуда, усл. ед");
xlabel("Номер отсчета");

figure(1);
plot(t*fs,imp_har1);
grid on;
ylabel("Амплитуда, усл. ед");
xlabel("Номер отсчета");