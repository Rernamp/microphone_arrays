clear;
close all;
%WideBand BeamForming
N = 16; 
J = 10;
N_f = 50;
fs = 8000;
phi_const = 40;
teta_const = 30;

c = 343;
p_loc = gen_place_el(4,4,0.04,0.04,1)';
w_nm = rand(1,N*J);


w_nm = ones(1,N*J);

[B,BB] = plot_bp_for_place(p_loc,phi_const,teta_const,N_f,w_nm,fs,J,N);
% WW = plot_wb_bp_freq_3d(w_nm, N, J, p_loc, N_f, phi_const, teta_const );