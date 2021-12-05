close all;
clear;

fs = 30000;
time = 2;
N = 16;
t = 0:1/fs:time;
f = 10;
y = sin(2*pi*f*t) + sin(2*pi*(f+10)*t);
maxi = max(y);
mini = min(y);
del_A = 2*(maxi - mini)/(2^N);

y_a = floor(y./del_A);
qe = 0;

for i = 1:length(y_a)
    d = y_a(i);
    if (d >= qe)
        pdm_y(i) = 1;
    end
    if (d < qe)
        pdm_y(i) = 0;
    end
    
    qe = pdm_y(i) - d + qe;
end

%%

figure()
plot(t,pdm_y)

%%
N = 16;
y_pcm = zeros(1, length(y_a)-15);
for i = 1:length(y_a)-15
    
    for k = 0:15
        y_pcm(i) = y_pcm(i) + pdm_y(i+k)*2^(k);
    end
    
    
    
end

%%
figure()
plot(t(1:length(y_a)-15),y_pcm)