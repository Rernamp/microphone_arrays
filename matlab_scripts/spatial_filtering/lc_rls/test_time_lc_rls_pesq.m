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
    sig_in_MR_noise(k,:) = y_shift+sig';
    sig_in_MR_sig(k,:) = sig;
end


time = 0:1/fs:(length(noise)-1)/fs;

L_max = 64;
L_min = 16;
L = L_min:2:L_max; 
i = 1;
%%
for l = L_min:2:L_max
    [y_cl,y,W_n] = spat_filt_wb_time_lc_rls_pesq(sig_in_MR_sig,sig_in_MR_noise, l, K);
    k_l(i,1:2) = pesqbin(y_cl,y,fs,'nb');
    i=i+1;
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
        sig_in_MR_noise(n,:) = y_shift+sig';
        sig_in_MR_sig(n,:) = sig;
        
    end
    [y_cl,y,W_n] = spat_filt_wb_time_lc_rls_pesq(sig_in_MR_sig,sig_in_MR_noise, L_k, p);
    k_K(p-k_min+1,1:2) = pesqbin(y_cl,y,fs,'nb');
end
%%
l = L_min:2:L_max;

figure()

plot(l,k_l(:,1))
grid on
title("PESQ для LC RLS")
xlabel("J , K = 4")
ylabel("PESQ")

K = k_min:k_max;

figure()

plot(K,k_K(:,1))
grid on
title("PESQ для LC RLS")
xlabel("K, J = 32")
ylabel("PESQ")

save('PESQ_RLS.mat','l','K','k_l','k_K');
