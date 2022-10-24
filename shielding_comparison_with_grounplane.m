% Comparison between ideal shield and discretized shield

%Data is taken from simulation files 
% "C:\Users\objoerkqvist\OneDrive - ETH Zurich\Documents\CST\Harvesting coil designs\vertical shielded coil v6\vertical_shielded_coil_v6_shield.cst"
% "C:\Users\objoerkqvist\OneDrive - ETH Zurich\Documents\CST\Harvesting coil designs\vertical shielded coil v6\vertical_shielded_coil_v6_GP_currents.cst"
close all

path_shield = "C:\Users\objoerkqvist\OneDrive - ETH Zurich\Documents\CST\Harvesting coil designs\vertical shielded coil v6\v6_shield_field_strength_cut_x=2.txt";
path_second_shield = "C:\Users\objoerkqvist\OneDrive - ETH Zurich\Documents\CST\Harvesting coil designs\vertical shielded coil v6\v6_second_shield_field_strength_cut_x=2.txt";
path_GP = "C:\Users\objoerkqvist\OneDrive - ETH Zurich\Documents\CST\Harvesting coil designs\vertical shielded coil v6\v6_GP_field_strength_cut_x=2.txt";

data_shield = importdata(path_shield);
data_second_shield = importdata(path_second_shield);
data_GP = importdata(path_GP);

Y_shield = reshape(data_shield.data(:,2),[numel(unique(data_shield.data(:,2))),numel(data_shield.data(:,2))/numel(unique(data_shield.data(:,2)))]);
Z_shield = reshape(data_shield.data(:,3),[numel(unique(data_shield.data(:,2))),numel(data_shield.data(:,2))/numel(unique(data_shield.data(:,2)))]);
Hx_shield = reshape(data_shield.data(:,4),[numel(unique(data_shield.data(:,2))),numel(data_shield.data(:,2))/numel(unique(data_shield.data(:,2)))]);
Hy_shield = reshape(data_shield.data(:,6),[numel(unique(data_shield.data(:,2))),numel(data_shield.data(:,2))/numel(unique(data_shield.data(:,2)))]);
Hz_shield = reshape(data_shield.data(:,8),[numel(unique(data_shield.data(:,2))),numel(data_shield.data(:,2))/numel(unique(data_shield.data(:,2)))]);
ABS_shield = sqrt(Hx_shield.^2 + Hy_shield.^2 + Hz_shield.^2)*0.6875;
ABS_shield_db = 20*log10(abs(ABS_shield));


Y_second_shield = reshape(data_second_shield.data(:,2),[numel(unique(data_second_shield.data(:,2))),numel(data_second_shield.data(:,2))/numel(unique(data_shield.data(:,2)))]);
Z_second_shield = reshape(data_second_shield.data(:,3),[numel(unique(data_second_shield.data(:,2))),numel(data_second_shield.data(:,2))/numel(unique(data_shield.data(:,2)))]);
Hx_second_shield = reshape(data_second_shield.data(:,4),[numel(unique(data_second_shield.data(:,2))),numel(data_second_shield.data(:,2))/numel(unique(data_shield.data(:,2)))]);
Hy_second_shield = reshape(data_second_shield.data(:,6),[numel(unique(data_second_shield.data(:,2))),numel(data_second_shield.data(:,2))/numel(unique(data_shield.data(:,2)))]);
Hz_second_shield = reshape(data_second_shield.data(:,8),[numel(unique(data_second_shield.data(:,2))),numel(data_second_shield.data(:,2))/numel(unique(data_shield.data(:,2)))]);
ABS_second_shield = sqrt(Hx_second_shield.^2 + Hy_second_shield.^2 + Hz_second_shield.^2)*0.667;
ABS_second_shield_db = 20*log10(abs(ABS_second_shield));


Y_GP = reshape(data_GP.data(:,2),[numel(unique(data_GP.data(:,2))),numel(data_GP.data(:,2))/numel(unique(data_GP.data(:,2)))]);
Z_GP = reshape(data_GP.data(:,3),[numel(unique(data_GP.data(:,2))),numel(data_GP.data(:,2))/numel(unique(data_GP.data(:,2)))]);
Hx_GP = reshape(data_GP.data(:,4),[numel(unique(data_GP.data(:,2))),numel(data_GP.data(:,2))/numel(unique(data_GP.data(:,2)))]);
Hy_GP = reshape(data_GP.data(:,6),[numel(unique(data_GP.data(:,2))),numel(data_GP.data(:,2))/numel(unique(data_GP.data(:,2)))]);
Hz_GP = reshape(data_GP.data(:,8),[numel(unique(data_GP.data(:,2))),numel(data_GP.data(:,2))/numel(unique(data_GP.data(:,2)))]);
ABS_GP = sqrt(Hx_GP.^2 + Hy_GP.^2 + Hz_GP.^2);
ABS_GP_db = 20*log10(abs(ABS_GP));

figure(1);
surf(Y_GP,Z_GP,ABS_GP_db);
shading interp

figure(2);
surf(Y_shield,Z_shield,ABS_shield_db);
shading interp

figure(3);
surf(Y_second_shield,Z_second_shield,ABS_second_shield_db);
shading interp



uiPlot(Z_GP,Z_shield,Z_second_shield,ABS_GP_db,ABS_shield_db,ABS_second_shield_db);