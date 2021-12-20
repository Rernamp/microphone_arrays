%% PDM generation
close all;
% clear all
% 4MHz sampling frequency as given in question
fs = 4e6;
ts = 1/fs;

% 50 kHz signal frequency as given in question
f50k = 50e3;
t50k = 1/f50k;

% Let's generate data for 20 cycles of 50kHz
t = [0 : ts : 20*t50k]';

% original signal : 50kHz modulated by a pulse.
os = 0.5 * sin(f50k*2*pi*t) .* (t >= 5*t50k & t<= 15*t50k);

% PDM generation as Given in Wikipedia
% https://en.wikipedia.org/wiki/Pulse-density_modulation#Algorithm
pdm = zeros(length(os),1);
qe = 0;
for ii = 1 : length(os)
    if(os(ii) >= qe)
        pdm(ii) = 1;
    else
        pdm(ii) = -1;
    end
    qe = pdm(ii) - os(ii) + qe;
end

% change all the -1 to 0 to match format in the question.
pdm(pdm < 0) = 0;

hold on 
plot(pdm)
plot(0.5+os)
xlim([400 500])


[numer, denom] = cheby2(4, 23, 300e3/2e6);

filtered = filtfilt(numer, denom, pdm);
%%
plot(filtered)
