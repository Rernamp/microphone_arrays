clear;
%задаю начальные параметры

N = 10;
f = 2000;
c = 343;
lamda = c/f;
d = lamda/2;
fi = pi/6;

ksi = pi/2;

q= 0:5:180;

fif = 2*pi*q/360;
F = zeros(1,length(fif));

for i = 1:N
    F = F + ((exp((1i*2*pi*d*(i-1)*sin(fif-ksi))/(lamda)))/((N)^0.5));
end

figure(1);
plot(q,abs(F));
