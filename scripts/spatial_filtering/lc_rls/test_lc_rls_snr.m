clear;
close all;

c = 343;
f_re = 2000;%ВЮЯРНРЮ ДКЪ ЙНРНПНИ ЯРПНХРЭЯЪ лп
lamda = c/f_re;%ДКХМЮ БНКМШ 
K = 4;%ЙНК-БН ДЮРВХЙНБ
f_up = 4000;    %БЕПУМЪЪ ВЮЯРНРЮ 
f_down = 0;  %МХФМЪЪ ВЮЯРНРЮ
f_den = 0.5;    %ЬЮЦ ВЮЯРНР
J = 30;
d = 0.04;%ПЮЯЯРНЪМХЕ ЛЕФДС ПЕЬЕРЙЮЛХ

thetaad = 0;       %ОПНОСЯЙЮЕЛШИ СЦНК
thetaan = 60;       %ГЮЦКСЬЮЕЛШИ СЦНК 

[noise,fs] = audioread('nois.wav');%ЯВХРШБЮЧ ЬСЛ

[sig,fs] = audioread('speech_dft_8kHz.wav');%ЯВХРШБЮЧ ОНКЕГМШИ ЯХЦМЮК

t = 2;%БПЕЛЪ ЯХЦМЮКЮ
%%
sig = sig(1:t*fs);
noise = noise(1:t*fs);

osh_in = mean(sig.^2)/mean(noise.^2);
sig_f = fft(sig);
sig_f(1) = 0;
sig = ifft(sig_f);

sig_in_MR_noise = zeros(K, length(sig));
sig_in_MR_sig = zeros(K, length(sig));

tau = (d*sind(thetaan))/(c);
%%
for k = 1:K
    [y_shift] = shift( noise, tau*(k-1), fs);
    sig_in_MR_noise(k,:) = y_shift;
    sig_in_MR_sig(k,:) = sig;
end


time = 0:1/fs:(length(noise)-1)/fs;

L_max = 64;
L_min = 16;
L = L_min:2:L_max; 
i = 1;
for l = L_min:2:L_max
    [y_noise,W_n] = spat_filt_wb_time_lc_rls(sig_in_MR_noise, l, K);
    [y_sig,W_s] = spat_filt_wb_time_lc_rls(sig_in_MR_sig, l, K);
    osh_l(i) = mean(y_sig.^2)/mean(y_noise.^2);
    i = i + 1;
end
%%
k_max = 8;
k_min = 2;
L_k = 32;
for p = k_min:k_max
    sig_in_MR_noise = zeros(p, length(sig));
    sig_in_MR_sig = zeros(p, length(sig));
    for n = 1:p
        [y_shift] = shift( noise, tau*(n-1), fs);
        sig_in_MR_noise(n,:) = y_shift;
        sig_in_MR_sig(n,:) = sig;
        
    end
    [y_noise,W_n] = spat_filt_wb_time_lc_rls(sig_in_MR_noise, L_k,p);
    [y_sig,W_s] = spat_filt_wb_time_lc_rls(sig_in_MR_sig, L_k, p);
    
    %y_noise = y_noise(41:end);
    %y_sig = y_sig(41:end);
    
    y_noise = y_noise;
    y_sig = y_sig;
    osh_k(p-k_min+1) = mean(y_sig.^2)/mean(y_noise.^2)
   
end
%%
l = L_min:2:L_max;

figure()
plot(l,db(osh_l-osh_in))
ylabel('ОСШ')
xlabel('Порядок фильтра, K = 4')
title('ОСШ взависимости от порядка фильтра')
grid on
%%
p = k_min:k_max;

figure()
plot(p,db(osh_k-osh_in))
ylabel('ОСШ')
xlabel('Число микрофонов, J = 32')
title('ОСШ взависимости от числа элементов МР')
grid on

