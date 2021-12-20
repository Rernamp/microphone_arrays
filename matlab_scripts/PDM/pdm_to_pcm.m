clear;
close all;

fileID = fopen('logdata.txt','r');
fs = 8e3;
A = fscanf(fileID,'%c');

k = length(A)*0.4;

data = A(1:k);

data_int = [];

for i = 1:length(data)/2
    data_int = [data_int hex2dec(data(2*(i)-1:2*i))];
end

data_bin = dec2bin(data_int);

data_pdm = double(data_bin) - 48;

data_1 = data_pdm(1:2:end-1);
data_2 = data_pdm(2:2:end);

x_ff = fft(data_1);
x_ff(1) = 0;
data_1 = ifft(x_ff);

[numer, denom] = cheby2(4, 23, 8e3/512e3);



pcm = filtfilt(numer, denom, data_1);
pcm = pcm(1:64:end);
figure()
plot(pcm)



