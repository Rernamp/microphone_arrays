clear; close all

load PESQ_LMS.mat
pesq_lms = k_l(:,1);
l_pesq_lms = l;
%%
load PESQ_LMS_move.mat
pesq_lms_move = k_l(:,1);
l_pesq_lms_move = l;
figure()
hold on
plot(l_pesq_lms(1:21),pesq_lms(1:21),'LineWidth',3)
plot(l_pesq_lms_move,pesq_lms_move,'LineWidth',3)
grid on
xlabel("Порядок фильтра")
ylabel("Оценка PESQ ")
legend("Исходная МР","МР со смещённым элементом",'Location','southeast')


load PESQ_RLS.mat
pesq_rls = k_l(:,1);
l_pesq_rls = l;
%%
load PESQ_RLS_move.mat
pesq_rls_move = k_l(:,1);
l_pesq_rls_move = l;
figure()
hold on
plot(l_pesq_rls,pesq_rls,'LineWidth',3)
plot(l_pesq_rls_move,pesq_rls_move,'LineWidth',3)
grid on
xlabel("Порядок фильтра")
ylabel("Оценка PESQ ")
legend("Исходная МР","МР со смещённым элементом",'Location','southeast')
%%
load SNR_RLS.mat
snr_rls = SNR_L(1:21);
l_snr_rls = l(1:21);

load SNR_RLS_move.mat
snr_rls_move = SNR_L;
l_snr_rls_move = l;
figure()
hold on
plot(l_snr_rls,snr_rls,'LineWidth',3)
plot(l_snr_rls_move,snr_rls_move,'LineWidth',3)
grid on
xlabel("Порядок фильтра")
ylabel("Выигрыш ОСШ, дБ")
legend("Исходная МР","МР со смещённым элементом",'Location','southeast')
%%
load SNR_LMS.mat
snr_LMS = SNR_L(1:21);
l_snr_LMS = l(1:21);

load SNR_LMS_move.mat
snr_LMS_move = SNR_L;
l_snr_LMS_move = l;
figure()
hold on
plot(l_snr_LMS,snr_LMS,'LineWidth',3)
plot(l_snr_LMS_move,snr_LMS_move,'LineWidth',3)
grid on
xlabel("Порядок фильтра")
ylabel("Выигрыш ОСШ, дБ")
legend("Исходная МР","МР со смещённым элементом",'Location','southeast')

%%

K = 9; 
d = 0.04;
p_loc = gen_place_el(sqrt(K),sqrt(K),d,d,1)';
p_cil = p_loc;
p_loc(3,6) = p_loc(3,6) + d/10; 
p_loc(2,6) = p_loc(2,6) - d/20; 

figure()
hold on
plot(p_loc(2,:),p_loc(3,:),'d','LineWidth',5, 'MarkerEdgeColor','r')
plot(p_cil(2,:),p_cil(3,:),'d','LineWidth',5, 'MarkerEdgeColor','b')
xline(0)
yline(0)
legend("МР со сдвинутым элементом","Прямоугольная равномерная МР")
xlim([-d*1.5 d*1.5])
ylim([-d*1.5 d*1.5])
grid on



