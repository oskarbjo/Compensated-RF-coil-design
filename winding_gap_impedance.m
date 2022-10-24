close all; clear all; clc

gaps = ["Winding gap = 0.1 mm","Winding gap = 0.48 mm","Winding gap = 0.86 mm","Winding gap = 1.24 mm","Winding gap = 1.62 mm","Winding gap = 2 mm"];

load winding_gap_impedance.mat;
freq = 1000*winding_gap_impedance(:,1).ParametersLinear_Spiral_cst_torrus_N10Linear_Spiral_cst_torrus_;
R1 = (winding_gap_impedance(:,2).VarName2);
N = numel(R1)/6;


figure(1);
semilogy(freq(1:N),R1(1:N),'linewidth',2)
hold on;
semilogy(freq(1:N),R1(N+1:2*N),'linewidth',2);
semilogy(freq(1:N),R1(2*N+1:3*N),'linewidth',2);
semilogy(freq(1:N),R1(3*N+1:4*N),'linewidth',2);
semilogy(freq(1:N),R1(4*N+1:5*N),'linewidth',2);
semilogy(freq(1:N),R1(5*N+1:6*N),'linewidth',2);
grid on;
xlim([0,400])
ylim([1e-1,1e6])
ylabel('Real impedance [Ohm]');
xlabel('Frequency [MHz]')
legend([gaps(:)])
