clear;
close all

load PESQ_LMS.mat
load SNR_LMS.mat

PESQ_LMS_K = k_K(:,1);
PESQ_LMS_L = k_l(:,1);

SNR_LMS_K = SNR_K;
SNR_LMS_L = SNR_L;

load PESQ_RLS.mat
load SNR_RLS.mat

PESQ_RLS_K = k_K(:,1);
PESQ_RLS_L = k_l(:,1);

SNR_RLS_K = SNR_K;
SNR_RLS_L = SNR_L;

L = l;
clear SNR_K SNR_L k_K k_l l;

figure()
hold on 
plot(L,PESQ_LMS_L)
plot(L,PESQ_RLS_L)
xlabel('Порядок фильтра')
ylabel('Оценка PESQ')
legend('LMS','RLS','Location','southeast')
grid on

figure()
hold on 
plot(L,SNR_LMS_L)
plot(L,SNR_RLS_L)
legend('LMS','RLS','Location','southeast')
xlabel('Порядок фильтра')
ylabel('Выигрыш в ОСШ, дБ')
grid on

figure()
hold on 
plot(K,SNR_LMS_K)
plot(K,SNR_RLS_K)
legend('LMS','RLS','Location','southeast')
xlabel('Число микрофонов')
ylabel('Выигрыш в ОСШ, дБ')
grid on

figure()
hold on 
plot(K,PESQ_LMS_K)
plot(K,PESQ_RLS_K)
legend('LMS','RLS','Location','southeast')
xlabel('Число микрофонов')
ylabel('Оценка PESQ')
grid on