clear;
clc;
close all;


[nois, Fs_n] = audioread('nois.wav');
[sign, Fs_s] = audioread('dft_voice_8kHz.wav');

t_duration = 2;%длительность сигнала
t = 0:1/Fs_s:t_duration-1/Fs_s;
NSampPerFrame = length(t);

nois(length(t)+1:end) = [];
sign(length(t)+1:end) = [];



%определяю задержку
n_tau = 2; %количество задерженных элементов
N_one = 1;%число отсчётов в одном бине




N_onezona = fix(length(sign)/(2*N_one))-1;      %кол-во целых диапазонов в первой зоне Найквиста
N_end1 = (fix(length(sign)/(2*N_one))*N_one)+2; %начало последнего диапазон
N_end2 = (length(sign)/2)+1;                    %конец последнего диапазона
N_ser = (N_end1 + N_end2)/2;                    % центральный отсчёт последнего диапазона
%для каждой решетки смешаю на своё число отсчётов

%создаю весовые коэффициенты для последующего превращения) в ИХ

i = 0:N_onezona;
f_tau_i = ((1+1+N_one-1)/2 + i*N_one);%центральная частота

    
w_one_zona = [exp(1i*2*pi*n_tau*f_tau_i/length(sign)) exp(1i*2*pi*n_tau*N_ser/length(sign))];
    
w_two_zona = conj(w_one_zona(end:-1:1));
w_w = [1 w_one_zona w_two_zona];

imp_har = ifft(w_w,length(w_w),2);


imp_har = fftshift(imp_har);


y_out = filter(imp_har,1,sign).';

%%
set(0,'DefaultAxesFontSize',14,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',14,'DefaultTextFontName','Times New Roman'); 


figure(1);
hold on;
grid on;
plot(1:length(imp_har),imp_har);
ylabel('Амплитуда, усл.ед');
xlabel('Время, с');

figure(2);
hold on;
grid on;
plot(t*Fs_n,y_out);
plot(t*Fs_n,sign);
legend("сдвиг","исход")
ylabel('Амплитуда, усл.ед');
xlabel('Время, с');

