clear;close all;


K = 9; 
c = 343;

fs = 8000;

f = 2000;
phi_sig = 0;
teta_sig = 0;
phi_noise = 50;
teta_noise = 0;
c = 343;
d = c/(2*f);
d = 0.04;
L = 32;
N_f = 50; 
mu = 0.01;
c_noise =35;
p_loc = gen_place_el(sqrt(K),sqrt(K),d,d,1)';

[noise,fs] = audioread('nois.wav');
[signal,fs] = audioread('speech_dft_8kHz.wav');
osh_in = mean(signal.^2)/mean(noise.^2);
t = 2;
time = 0:1/fs:t-1/fs;
signal = signal(1:t*fs);
noise = noise(1:t*fs);

signal_shift = shift_plane(signal,phi_sig,teta_sig,p_loc,fs);
noise_shift = shift_plane(noise,phi_noise,teta_noise,p_loc,fs);


mu = 1;

time = 0:1/fs:(length(noise)-1)/fs;

L_max = 64;
L_min = 16;
L = L_min:L_max; 
i = 1;

for l = L_min:2:L_max
    [y_noise,W_n] = spat_filt_wb_time_lc_lms(noise_shift, l, K, mu);
    [y_sig,W_s] = spat_filt_wb_time_lc_lms(signal_shift, l, K, mu);

    osh_l(i) = mean(y_sig.^2)/mean(y_noise.^2);
    i = i+1;
   
   
end
%%
mu = 0.1;
k = [4 9 16 25];
L_k = 32;
for p = 1:length(k)
    
    
    p_loc = gen_place_el(sqrt(k(p)),sqrt(k(p)),d,d,1)'; 
    signal_shift = shift_plane(signal,phi_sig,teta_sig,p_loc,fs);
    noise_shift = shift_plane(noise,phi_noise,teta_noise,p_loc,fs);
        
    
    [y_noise,W_n] = spat_filt_wb_time_lc_lms(noise_shift, L_k,k(p), mu);
    [y_sig,W_s] = spat_filt_wb_time_lc_lms(signal_shift, L_k, k(p), mu);
    osh_k(p) = mean(y_sig.^2)/mean(y_noise.^2);
    
end
%%
l = L_min:2:L_max;
SNR_L = db(osh_l-osh_in);
figure()
plot(l,SNR_L)
ylabel('ОСШ')
xlabel('Порядок фильтра, K = 4')
title('ОСШ взависимости от порядка фильтра для LC RLS')
grid on
%%

SNR_K = db(osh_k-osh_in);
figure()
plot(k,SNR_K)
ylabel('ОСШ')
xlabel('Число микрофонов, J = 32')
title('ОСШ взависимости от числа элементов МР для LC RLS')
grid on

save('SNR_LMS.mat','l','K','SNR_K','SNR_L');
