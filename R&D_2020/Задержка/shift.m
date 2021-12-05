clear;
%задаю начальные параметры

N = 10;%кол-во решеток
c = 343;%скорость звука
f = 2000;%частота
lamda = c/f;
d = lamda/2;%рассто€ние между решетками
fi = 40;%угол,под который пришел шум

[noise,fs] = audioread('nois.wav');%считываю шум

t_duration = 2;%длительность сигнала
t = 0:1/fs:t_duration-1/fs;
NSampPerFrame = length(t);

hap = dsp.AudioPlayer('SampleRate',fs);
%.................................—игнал......................................
dftFileReader = dsp.AudioFileReader('dft_voice_8kHz.wav',...
    'SamplesPerFrame',NSampPerFrame);
sig = step(dftFileReader);

nois = noise(1:length(sig));%укорачиваю кол-во отсчЄтов шума,что бы совпадали

%определ€ю задержку
tau = (d*sind(fi))/(c);

dt = 1/fs;
n_tau = 300; %перевожу в отсчЄты

N_one = 50;%число отсчЄтов в одной последовательности


spec = fft(nois);%Ѕѕ‘ шума,ибо сигнал приходит при 0 угле и не притерпевает задержки ни на одной решетке

%создаю массивы,которые дальше буду испольховать

y_s_1 = [];% перва€ зона Ќайквиста
y_s_2 = [];%втора€ зона Ќайквиста
y_s_v = [];%сигнал на каждой решетке
N_onezona =fix(length(sig)/(2*N_one))-1%кол-во целых диапазонов в первой зоне Ќайквиста
N_end1 = (fix(length(sig)/(2*N_one))*N_one)+2;%начало последнего диапазон
N_end2 = (length(sig)/2)+1;%конец последнего диапазона
N_ser = (N_end1 + N_end2)/2; % центральный отсчЄт последнего диапазона
%дл€ каждой решетки смешаю на своЄ число отсчЄтов
for m = 1:N %цикл дл€ каждой решетки
    n_tau_i = n_tau*(m-1);%кол-во задержаных отсчЄтов дл€ каждой рещетки
    for i = 0:N_onezona% цикл,который производит задержку шума
        
        
        f_tau_i = ((2+2+N_one-1)/2 + i*N_one);%отсчЄт,соответствующий центральной частоте дл€ каждого диапазона
        k1 = spec(2+i*N_one:1+N_one+i*N_one) .* exp(1i*2*pi*n_tau_i*f_tau_i/length(spec));%задержка
        y_s_1 = cat(1,y_s_1,k1);%присоедин€ю к пустому массиву
        
    end 

    y_s_1 = [y_s_1 (spec(N_end1:N_end2).*exp(1i*2*pi*n_tau_i*N_ser/length(spec)))];%добавл€ю отсчЄты,не вошедшие в основной массив
    y_s_1(end) = real(y_s_1(end));%так как отсчЄтов четное кол-во,то делаю реальным центральный отсчЄт 
    y_s_2(1:N_end2 - 2) = conj(y_s_1((length(sig)/2)-1:-1:(length(sig)/2)-N_end2 + 2));%комплексно сопр€гаю и отражаю отсчЄты дл€ второй зоны Ќайквиста
    y_f = [spec(1) y_s_1.' y_s_2];%соедин€ю в отсчЄты
    y_fdel = ifft(y_f); 
    y_fdel = y_fdel + sig.';%учитываю сигнал
        
    y_s_v = cat(1,y_s_v,y_fdel);%присоедин€ю полученный сигнал и так дл€ каждой решетки
    
    %обнул€ю массивы
    
    y_s_1 = [];
    y_s_2 = [];
    y_f = [];
    y_fdel = [];
    
end
 

%%
set(0,'DefaultAxesFontSize',14,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',14,'DefaultTextFontName','Times New Roman'); 

k = 100;
N_zad = 3;
nois_T_zad = [nois(N_zad+1:end,1); zeros(N_zad,1)];
sig_T_zad = nois + sig;
figure(1);
hold on;
grid on;
plot(t,sig_T_zad);
plot(t,y_s_v(2,:));

ylabel('јмплитуда, усл.ед');
xlabel('¬рем€, с');