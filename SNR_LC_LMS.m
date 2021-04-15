clear;close all;

c = 343;
f_re = 4000;
lamda = c/f_re;
K = 4;
f_up = 4000;   
f_down = 0;  
f_den = 0.5;   
L = 30;
d = 0.04;

thetaad = 0;     
thetaan = 60;      

[noise,fs] = audioread('nois.wav');

[sig,fs] = audioread('speech_dft_8kHz.wav');

t = 2;
%%
sig = sig(1:t*fs);
noise = noise(1:t*fs);

sig_in_MR_noise = zeros(K, length(sig));
sig_in_MR_sig = zeros(K, length(sig));

tau = (d*sind(thetaan))/(c);
osh_in = mean(sig.^2)/mean(noise.^2);
%%
for k = 1:K
    [y_shift] = shift( noise, tau*(k-1), fs);
    sig_in_MR_noise(k,:) = y_shift;
    sig_in_MR_sig(k,:) = sig;
end

mu = 0.05;

time = 0:1/fs:(length(noise)-1)/fs;

L_max = 64;
L_min = 16;
L = L_min:L_max; 


for l = L_min:L_max
    [y_noise,W_n] = func_Frost(sig_in_MR_noise, l, K, mu);
    [y_sig,W_s] = func_Frost(sig_in_MR_sig, l, K, mu);
    osh_l(l-L_min+1) = mean(y_sig.^2)/mean(y_noise.^2)
   
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
    [y_noise,W_n] = func_Frost(sig_in_MR_noise, L_k,p, mu);
    [y_sig,W_s] = func_Frost(sig_in_MR_sig, L_k, p, mu);
    osh_k(p-k_min+1) = mean(y_sig.^2)/mean(y_noise.^2);
    
end
%%
l = L_min:L_max;

figure()

plot(l,db(osh_l-osh_in))
grid on
title("ОСШ от порядка фильтра")
xlabel("J, K = 4")
ylabel("ОСШ")

p = k_min:k_max;

figure()

plot(p,db(osh_k-osh_in))
grid on
title("ОСШ от числа элементов")
xlabel("K, J = 32")
ylabel("ОСШ")