clear;
%задаю начальные параметры

set(0,'DefaultAxesFontSize',14,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',14,'DefaultTextFontName','Times New Roman'); 

N = 10;%кол-во решеток
c = 343;%скорость звука
f = 2000;%частота
lamda = c/f;
d = lamda/2;%рассто€ние между решетками
fi = 40;%угол,под который пришел шум
N_one = 20;%число отсчЄтов в одной последовательности

[noise,fs] = audioread('nois.wav');%считываю шум

t_duration = 2;%длительность сигнала
t = 0:1/fs:t_duration-1/fs;
NSampPerFrame = length(t);

hap = dsp.AudioPlayer('SampleRate',fs);
%.................................—игнал......................................
dftFileReader = dsp.AudioFileReader('dft_voice_8kHz.wav',...
    'SamplesPerFrame',NSampPerFrame);
sig = step(dftFileReader);
f = (0:length(sig(:,1))-1)/t_duration;

nois = noise(1:length(sig));%укорачиваю кол-во отсчЄтов шума,что бы совпадали

spec_nois = fft(nois);
spec_sig = fft(sig);

ahc_nois = mag2db(abs(spec_nois));
ahc_sig = mag2db(abs(spec_sig));

% figure('Units', 'normalized', 'OuterPosition', [0 0 1 1]);
%%
figure(1);
plot(f,ahc_sig);
grid on;
xlabel('„астота, Hz');
ylabel("BP, дЅ")

figure(2);
plot(f,ahc_nois);
grid on;

xlabel('„астота, Hz');
ylabel("BP, дЅ")