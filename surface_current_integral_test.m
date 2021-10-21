%% Surface current integral test
clc;
close all;

path="C:\Users\objoerkqvist\Documents\CST\Harvesting coil designs\surface current integral test\exported fields\top_plane_real.txt";
data = importdata(path);
data=struct2cell(data);
x_ind = 1;
y_ind = 2;
z_ind = 3;
Hx_ind = 4;
Hy_ind = 6;
Hz_ind = 8;


x = data{1}(:,x_ind);
y = data{1}(:,y_ind);
z = data{1}(:,z_ind);
Hx = data{1}(:,Hx_ind);
Hy = data{1}(:,Hy_ind);
Hz = data{1}(:,Hz_ind);

zeroVector = zeros(length(Hx),1);
oneVector = ones(length(Hx),1);
J = cross([Hx,Hy,Hz],[zeroVector,zeroVector,oneVector]);
Jx = J(:,1);
Jy = J(:,2);
Jz = J(:,3);
Jabs = sqrt(Jx.^2 + Jy.^2 + Jz.^2);

%calculate x y grid, assuming perfect square grid
[dummie,xmin_ind] = min(x);
[dummie,xmax_ind] = max(x);
[X,Y] = meshgrid(x(xmin_ind:xmax_ind),x(xmin_ind:xmax_ind));
JX = reshape(Jx,size(X))';
JY = reshape(Jy,size(X))';
J_Abs = sqrt(JX.^2 + JY.^2);

figure();
quiver(X,Y,JX,JY);


%% Calculate integral

%Find x = 0 index
zero_ind=find(X==0)/length(X);
zero_ind = zero_ind(end);

%Find delta r
dr=0.001*(X(1,2)-X(1,1));

%Integrate JX along x=0
current_integral = cumsum(JX(:,zero_ind))*dr;

figure();
plot(JX(:,zero_ind));

%Cumsum should reach ~1/2 of total current of 1A (the rest of the current 
%is confined to bottom half of conductor, assuming trace width >> thickness). 
%Seems to add up!
figure();
plot(current_integral);

