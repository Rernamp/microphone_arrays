%snr
clear; close all

load SNR_LMS.mat
SNR_LMS_K = SNR_K;
SNR_LMS_L = SNR_L;
load SNR_RLS.mat
SNR_RLS_K = SNR_K;
SNR_RLS_L = SNR_L;
k = [4 9 16 25];

%snr for K - number el
figure()
hold on
plot(k,SNR_LMS_K,'LineWidth',3)
plot(k,SNR_RLS_K,'LineWidth',3)
grid on
xlabel("Число микрофонов")
ylabel("Выигрыш ОСШ, дБ")
legend("LC NLMS","LC RLS",'Location','northwest')

%snr for L - por filt
figure()
hold on
plot(l,SNR_LMS_L,'LineWidth',3)
plot(l,SNR_RLS_L,'LineWidth',3)
grid on
ylabel("Выигрыш ОСШ, дБ")
xlabel("Порядок фильтра")
legend("LC NLMS","LC RLS",'Location','northwest')

clear
load PESQ_LMS.mat
l_lms = l;
PESQ_LMS_K = k_K;
PESQ_LMS_L = k_l;
load PESQ_RLS.mat
PESQ_RLS_K = k_K;
PESQ_RLS_L = k_l;
k = [4 9 16 25];

figure()
hold on
plot(k,PESQ_LMS_K(:,1),'LineWidth',3)
plot(k,PESQ_RLS_K(:,1),'LineWidth',3)
grid on
xlabel("Число микрофонов")
ylabel("Оценка PESQ")
legend("LC NLMS","LC RLS",'Location','northwest')

%snr for L - por filt
figure()
hold on
plot(l_lms(1:21),PESQ_LMS_L(1:21,1),'LineWidth',3)
plot(l,PESQ_RLS_L(:,1),'LineWidth',3)
grid on
ylabel("Оценка PESQ")
xlabel("Порядок фильтра")
legend("LC NLMS","LC RLS",'Location','northwest')


