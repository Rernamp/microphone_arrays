clear; close all
load PESQ_LMS.mat
PESQ_LMS_phi = k_phi;
PESQ_LMS_teta = k_teta;
phi_lms = phi;
teta_lms = teta;


load PESQ_RLS.mat
PESQ_RLS_phi = k_phi;
PESQ_RLS_teta = k_teta;
phi_rls = phi;
teta_rls = teta;

figure()
hold on
plot(phi_lms,PESQ_LMS_phi(:,1),'LineWidth',3)
plot(phi_rls,PESQ_RLS_phi(:,1),'LineWidth',3)

grid on;
ylabel("Оценка PESQ")
xlabel("Угол азимута, град")
legend("LC NLMS", "LC RLS")

figure()
hold on
plot(teta_lms,PESQ_LMS_teta(:,1),'LineWidth',3)
plot(teta_rls,PESQ_RLS_teta(:,1),'LineWidth',3)
grid on;
ylabel("Оценка PESQ")
xlabel("Угол подъёма, град")
legend("LC NLMS", "LC RLS")
%%
clear; 
load SNR_LMS.mat
SNR_LMS_phi = SNR_phi;
SNR_LMS_teta = SNR_teta;
phi_lms = phi;
teta_lms = teta;


load SNR_RLS.mat
SNR_RLS_phi = SNR_phi;
SNR_RLS_teta = SNR_teta;
phi_rls = phi;
teta_rls = teta;

figure()
hold on
plot(phi_lms,SNR_LMS_phi,'LineWidth',3)
plot(phi_rls,SNR_RLS_phi,'LineWidth',3)


grid on;
ylabel("Выигрыш ОСШ")
xlabel("Угол азимута, град")
legend("LC NLMS", "LC RLS")

figure()
hold on
plot(teta_lms,SNR_LMS_teta,'LineWidth',3)
plot(teta_rls,SNR_RLS_teta,'LineWidth',3)
grid on;
ylabel("Выигрыш ОСШ")
xlabel("Угол подъёма, град")
legend("LC NLMS", "LC RLS")