%% S11 plot

load multiturn_compensation_coils_S11_RI.mat;

S11_compl = S11_data(:,2)+S11_data(:,3)*i;
Z0 = 1;
Zin = Z0 * (1+S11_compl)./(1-S11_compl);

figure()
subplot(2,1,1)
plot(S11_data(:,1),real(Zin));
subplot(2,1,2);
plot(S11_data(:,1),imag(Zin));