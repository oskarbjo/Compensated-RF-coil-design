% RF harvester v4 VNA coupling measurements

clear all;
close all; clc;

PASSIVE_SHIELD_V4_1NF_TOP=sparameters("C:\Users\objoerkqvist\OneDrive - ETH Zurich\Documents\Measurements\RF harvester\Coupling top bottom v4\PASSIVE_SHIELD_V4_1NF_TOP.S1P");
PASSIVE_SHIELD_V4_1NF_BOTTOM=sparameters("C:\Users\objoerkqvist\OneDrive - ETH Zurich\Documents\Measurements\RF harvester\Coupling top bottom v4\PASSIVE_SHIELD_V4_1NF_BOTTOM.S1P");
PASSIVE_SHIELD_V4_56PF_TOP=sparameters("C:\Users\objoerkqvist\OneDrive - ETH Zurich\Documents\Measurements\RF harvester\Coupling top bottom v4\PASSIVE_SHIELD_V4_56PF_TOP.S1P");
PASSIVE_SHIELD_V4_56PF_BOTTOM=sparameters("C:\Users\objoerkqvist\OneDrive - ETH Zurich\Documents\Measurements\RF harvester\Coupling top bottom v4\PASSIVE_SHIELD_V4_56PF_BOTTOM.S1P");
PASSIVE_SHIELD_V4_68PF_TOP=sparameters("C:\Users\objoerkqvist\OneDrive - ETH Zurich\Documents\Measurements\RF harvester\Coupling top bottom v4\PASSIVE_SHIELD_V4_68PF_TOP.S1P");
PASSIVE_SHIELD_V4_68PF_BOTTOM=sparameters("C:\Users\objoerkqvist\OneDrive - ETH Zurich\Documents\Measurements\RF harvester\Coupling top bottom v4\PASSIVE_SHIELD_V4_68PF_BOTTOM.S1P");

S11_PASSIVE_SHIELD_V4_1NF_TOP(:)=PASSIVE_SHIELD_V4_1NF_TOP.Parameters(1,1,:);
S11_PASSIVE_SHIELD_V4_1NF_BOTTOM(:)=PASSIVE_SHIELD_V4_1NF_BOTTOM.Parameters(1,1,:);
figure()
plot(PASSIVE_SHIELD_V4_1NF_TOP.Frequencies,abs(S11_PASSIVE_SHIELD_V4_1NF_TOP));
hold on;
plot(PASSIVE_SHIELD_V4_1NF_BOTTOM.Frequencies,abs(S11_PASSIVE_SHIELD_V4_1NF_BOTTOM));
grid on;
xlim([100e6,150e6]);
ylim([0.4,1.2]);
legend('Unshielded side','Shielded side');
ylabel('S11 [linear]');
xlabel('Frequency');
title('1 nF shield capacitors');

S11_PASSIVE_SHIELD_V4_68PF_TOP(:)=PASSIVE_SHIELD_V4_68PF_TOP.Parameters(1,1,:);
S11_PASSIVE_SHIELD_V4_68PF_BOTTOM(:)=PASSIVE_SHIELD_V4_68PF_BOTTOM.Parameters(1,1,:);
figure()
plot(PASSIVE_SHIELD_V4_68PF_TOP.Frequencies,abs(S11_PASSIVE_SHIELD_V4_68PF_TOP));
hold on;
plot(PASSIVE_SHIELD_V4_68PF_BOTTOM.Frequencies,abs(S11_PASSIVE_SHIELD_V4_68PF_BOTTOM));
grid on;
xlim([100e6,150e6]);
ylim([0.4,1.2]);
legend('Unshielded side','Shielded side');
ylabel('S11 [linear]');
xlabel('Frequency');
title('68 nF shield capacitors');

S11_PASSIVE_SHIELD_V4_56PF_TOP(:)=PASSIVE_SHIELD_V4_56PF_TOP.Parameters(1,1,:);
S11_PASSIVE_SHIELD_V4_56PF_BOTTOM(:)=PASSIVE_SHIELD_V4_56PF_BOTTOM.Parameters(1,1,:);
figure()
plot(PASSIVE_SHIELD_V4_56PF_TOP.Frequencies,abs(S11_PASSIVE_SHIELD_V4_56PF_TOP));
hold on;
plot(PASSIVE_SHIELD_V4_56PF_BOTTOM.Frequencies,abs(S11_PASSIVE_SHIELD_V4_56PF_BOTTOM));
grid on;
xlim([100e6,150e6]);
ylim([0.4,1.2]);
legend('Unshielded side','Shielded side');
ylabel('S11 [linear]');
xlabel('Frequency');
title('56 nF shield capacitors');
