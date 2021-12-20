clear;close all;

c = 343;
f_re = 4000;
lamda = c/f_re;
K = 2;
d = 0.04;


thetaan = 60; 

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
    awgn(y_shift, 60);
    noise_shift(k,:) = y_shift;
end


t = 0:1/fs:2-1/fs;

a = [0.5 0.5];

X_fft = [fft(noise_shift(1,:)); fft(noise_shift(2,:))];

shift_fft = X_fft(1,:) - X_fft(2,:);

y_del = noise_shift(1,:) - noise_shift(2,:);
y_del = y_del ./noise_shift(1,:);

hold on
plot(y_del)

