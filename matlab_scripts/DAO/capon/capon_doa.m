clear;close all;

c = 343;
f_re = 4000;
lamda = c/f_re;
K = 2;
d = 0.04;


thetaan = 15; 

% [noise,fs] = audioread('sound_44.1kHz.wav');
[noise,fs] = audioread('sound.wav');

t = 2;

noise = noise(1:t*fs);



tau = (d*sind(thetaan))/(c)
% tau = 0.8/fs;

noise_shift = zeros(K, length(noise));

%%
for k = 1:K
    [y_shift] = shift( noise, tau*(k-1), fs);
    noise_shift(k,:) = y_shift;
end


t = 0:1/fs:2-1/fs;

a = [0.5 0.5];

X_fft = [fft(noise_shift(1,:)); fft(noise_shift(2,:))];

phi_find = [-90:90];

