load coil10mmradius1to10windingsS11.mat
load harvestcoil.mat;
load harvesting_coil_v7_param_sweep_R.mat;

f = freqGHz*1e9;
A = 0.01^2 * pi;
N = 1:10;
Z = (s11 + 1)./(1 - s11);
B1 = 1e-6;


emf = 2*pi*f*A*N*B1;


% freqGHz=coil10mmradius1to10windingss11data(1:1001,1);
% s11=   [coil10mmradius1to10windingss11data(1008-3:1008+1000-3,2)+1j*coil10mmradius1to10windingss11data(1008-3:1008+1000-3,3), ...
%         coil10mmradius1to10windingss11data(2013-3:2013+1000-3,2)+1j*coil10mmradius1to10windingss11data(2013-3:2013+1000-3,3), ...
%         coil10mmradius1to10windingss11data(3018-3:3018+1000-3,2)+1j*coil10mmradius1to10windingss11data(3018-3:3018+1000-3,3), ...
%         coil10mmradius1to10windingss11data(4023-3:4023+1000-3,2)+1j*coil10mmradius1to10windingss11data(4023-3:4023+1000-3,3), ...
%         coil10mmradius1to10windingss11data(1:1001,2)+1j*coil10mmradius1to10windingss11data(1:1001,3), ...
%         coil10mmradius1to10windingss11data(5027-3:5027+1000-3,2)+1j*coil10mmradius1to10windingss11data(5027-3:5027+1000-3,3), ...
%         coil10mmradius1to10windingss11data(6032-3:6032+1000-3,2)+1j*coil10mmradius1to10windingss11data(6032-3:6032+1000-3,3), ...
%         coil10mmradius1to10windingss11data(7036-3:7036+1000-3,2)+1j*coil10mmradius1to10windingss11data(7036-3:7036+1000-3,3), ...
%         coil10mmradius1to10windingss11data(8040-3:8040+1000-3,2)+1j*coil10mmradius1to10windingss11data(8040-3:8040+1000-3,3), ...
%         coil10mmradius1to10windingss11data(9044-3:9044+1000-3,2)+1j*coil10mmradius1to10windingss11data(9044-3:9044+1000-3,3)];

figure;
for(i=1:10)
    plot(f,(emf(:,i).^2)./(4*real(Z(:,i))))
    ylim([0,0.3]);
    hold on;
    % plot(f,R/1e8);
end
legend('1','2','3','4','5','6','7','8','9','10');


figure;
for(i=1:6)
%     plot((emf(:,i).^2)./(4*real(Z(:,i))))
    plot(f,real(Z(:,i)))
%     ylim([0,0.3]);
    hold on;
    % plot(f,R/1e8);
end
legend('1','2','3','4','5','6','7','8','9','10');



ff = harvestingcoilv7impedance(:,1)*1e9;
AA = 0.029*0.022;
N=3;
B1=1e-6;
emf = AA*N*ff*2*pi*B1;
RR = harvestingcoilv7impedance(:,2);
figure;
plot(ff/1e9,RR,'linewidth',1);
ylabel('Coil ESR [Ohm]');
xlabel('Frequency [MHz]');
grid on;

figure;
plot(ff,emf.^2 ./ RR,'linewidth',1)
ylabel('Power (EMF^2/ESR) [a.u.]');
xlabel('Frequency [MHz]');
grid on;




freq = freq * 1e9;
figure;
L = 0.029;
H = 0.014;
AA = H*L;
N=3;
B1=1e-6;
emf = AA*N*freq*2*pi*B1;
plot(freq,emf.^2 ./ R_N3_H14); hold on;

H = 0.016;
AA = H*L;
N=3;
emf = AA*N*freq*2*pi*B1;
plot(freq,emf.^2 ./ R_N3_H16);

H = 0.018;
AA = H*L;
N=3;
emf = AA*N*freq*2*pi*B1;
plot(freq,emf.^2 ./ R_N3_H18);

H = 0.020;
AA = H*L;
N=3;
emf = AA*N*freq*2*pi*B1;
plot(freq,emf.^2 ./ R_N3_H20);

H = 0.022;
AA = H*L;
N=3;
emf = AA*N*freq*2*pi*B1;
plot(freq,emf.^2 ./ R_N3_H22);
title('3 turns')
legend('Height = 14 mm','Height = 16 mm','Height = 18 mm','Height = 20 mm','Height = 22 mm');






figure;
L = 0.029;
H = 0.014;
AA = H*L;
N=2;
B1=1e-6;
emf = AA*N*freq*2*pi*B1;
plot(freq,emf.^2 ./ R_N2_H14); hold on;

H = 0.016;
AA = H*L;
N=2;
emf = AA*N*freq*2*pi*B1;
plot(freq,emf.^2 ./ R_N2_H16);

H = 0.018;
AA = H*L;
N=2;
emf = AA*N*freq*2*pi*B1;
plot(freq,emf.^2 ./ R_N2_H18);

H = 0.020;
AA = H*L;
N=2;
emf = AA*N*freq*2*pi*B1;
plot(freq,emf.^2 ./ R_N2_H20);


H = 0.022;
AA = H*L;
N=2;
emf = AA*N*freq*2*pi*B1;
plot(freq,emf.^2 ./ R_N2_H22);
title('2 turns')
legend('Height = 14 mm','Height = 16 mm','Height = 18 mm','Height = 20 mm','Height = 22 mm');





figure;
L = 0.029;
H = 0.014;
AA = H*L;
N=4;
B1=1e-6;
emf = AA*N*freq*2*pi*B1;
plot(freq,emf.^2 ./ R_N4_H14); hold on;

H = 0.016;
AA = H*L;
N=4;
emf = AA*N*freq*2*pi*B1;
plot(freq,emf.^2 ./ R_N4_H16);

H = 0.018;
AA = H*L;
N=4;
emf = AA*N*freq*2*pi*B1;
plot(freq,emf.^2 ./ R_N4_H18);

H = 0.020;
AA = H*L;
N=4;
emf = AA*N*freq*2*pi*B1;
plot(freq,emf.^2 ./ R_N4_H20);


H = 0.022;
AA = H*L;
N=4;
emf = AA*N*freq*2*pi*B1;
% plot(freq,emf.^2 ./ R_N4_H22);
title('4 turns')
legend('Height = 14 mm','Height = 16 mm','Height = 18 mm','Height = 20 mm','Height = 22 mm');




figure;
L = 0.029;
H = 0.014;
AA = H*L;
N=5;
B1=1e-6;
emf = AA*N*freq*2*pi*B1;
plot(freq,emf.^2 ./ R_N5_H14); hold on;

H = 0.016;
AA = H*L;
N=5;
emf = AA*N*freq*2*pi*B1;
plot(freq,emf.^2 ./ R_N5_H16);

H = 0.018;
AA = H*L;
N=5;
emf = AA*N*freq*2*pi*B1;
plot(freq,emf.^2 ./ R_N5_H18);

H = 0.020;
AA = H*L;
N=5;
emf = AA*N*freq*2*pi*B1;
plot(freq,emf.^2 ./ R_N5_H20);


H = 0.022;
AA = H*L;
N=5;
emf = AA*N*freq*2*pi*B1;
plot(freq,emf.^2 ./ R_N5_H22);
title('5 turns');
legend('Height = 14 mm','Height = 16 mm','Height = 18 mm','Height = 20 mm','Height = 22 mm');



figure;
H = 0.022;
AA = H*L;
N=2;
emf = AA*N*freq*2*pi*B1;
plot(freq,emf.^2 ./ R_N2_H22);
hold on;

AA = H*L;
N=3;
emf = AA*N*freq*2*pi*B1;
plot(freq,emf.^2 ./ R_N3_H22);

AA = H*L;
N=4;
emf = AA*N*freq*2*pi*B1;
% plot(freq,emf.^2 ./ R_N4_H22);

AA = H*L;
N=5;
emf = AA*N*freq*2*pi*B1;
plot(freq,emf.^2 ./ R_N5_H22);

title('H = 22, N sweep');
legend('2 Turns','3 Turns','5 Turns')
