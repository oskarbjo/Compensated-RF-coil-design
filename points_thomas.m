% Prep data for Thomas S


load streamlines.mat;
L1 = 15;
wg = 0.25;
spread_factor = 1.2;

line1_x = (-L1/2:0.1:L1/2)';
line1_y = zeros(size(line1_x));
line1_z = zeros(size(line1_x));

phi = (-pi/2:0.02:pi/2)';
R = 5/2;
line2_x = line1_x(end) + R * cos(phi);
line2_y = zeros(size(line2_x));
line2_z = R + R * sin(phi);

line3_x = fliplr(line1_x')';
line3_y = fliplr(line1_y);
line3_z = line2_z(end)*ones(size(line1_z));

line4_x = line3_x(end) - R*cos(phi);
line4_y = (phi+pi/2)/pi * wg * spread_factor;
line4_z = R - R * sin(phi);

output=[line1_x, line1_y, line1_z];
output=[output;[line2_x, line2_y, line2_z]];
output=[output;[line3_x, line3_y, line3_z]];
output=[output;[line4_x, line4_y, line4_z]];

N = length(output);
output1 = output;
for i = 1:8
   
   output_translated = output1;
   output_translated(:,2) = output_translated(:,2) + i*output1(N,2);
   output = [output; output_translated]; 
    
end

figure();
for i = 1:length(streamLineList)
    streamLineList{i} = [streamLineList{i},zeros(length(streamLineList{i}),1)];
    plot3(streamLineList{i}(:,1),streamLineList{i}(:,2),streamLineList{i}(:,3));
    hold on;
end

plot3(output(:,1),output(:,2),output(:,3));

%%

load Thomas_Data.mat

figure();
for i = 1:length(streamLineList)
    plot3(streamLineList{i}(:,1),streamLineList{i}(:,2),streamLineList{i}(:,3));
    hold on;
    axis equal;
end

plot3(output(:,1),output(:,2),output(:,3));
