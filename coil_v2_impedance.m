load vertical_shielded_coil_full_simulation_impedance_vs_freq.txt
data = vertical_shielded_coil_full_simulation_impedance_vs_freq;



figure();
plot(data(:,1)*1000,data(:,2));
hold on;
plot(data(:,1)*1000,data(:,3));
grid on;
legend('Real part','Imaginary part');
xlabel('Frequency [MHz]');
ylabel('Impedance [Ohm]');

L = data(241,3)/(data(241,1)*2*pi*1e9)w*