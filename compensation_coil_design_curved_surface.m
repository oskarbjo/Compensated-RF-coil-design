%% Compensation coil design curved surface

path="C:\Users\objoerkqvist\Documents\CST\Harvesting coil designs\warped shielded coil\Hfield_curved_GP_2.txt";
data = importdata(path);
data=struct2cell(data);

%%

%Curve shape: -0.01*x^2

x = linspace(min(data{1}(:,1)),max(data{1}(:,1)),length(unique(data{1}(:,1))));
y = linspace(min(data{1}(:,2)),max(data{1}(:,2)),length(unique(data{1}(:,2))));
z = linspace(min(data{1}(:,3)),max(data{1}(:,3)),length(unique(data{1}(:,3))));

[X,Y,Z]= meshgrid(x,y,z);
X=permute(X,[2,1,3]);
Y=permute(Y,[2,1,3]);


HX = reshape(data{1}(:,4),size(X));
% HX = permute(HX,[2,1,3]);
HY = reshape(data{1}(:,6),size(Y));
% HY = permute(HY,[2,1,3]);
HZ = reshape(data{1}(:,8),size(Z));

figure();
zplane = 10;
quiver(X(:,:,zplane),Y(:,:,zplane),HX(:,:,zplane),HY(:,:,zplane));



%% Find vectors that lie on surface z=-0.01*x^2 (stepping from underneath the surface, field should be zero until surface is reached)

Hx_surf = zeros(length(unique(data{1}(:,1))),length(unique(data{1}(:,2))));
Hy_surf = zeros(length(unique(data{1}(:,1))),length(unique(data{1}(:,2))));
Hz_surf = zeros(length(unique(data{1}(:,1))),length(unique(data{1}(:,2))));

figure();

for z_level = length(unique(data{1}(:,3))):-1:1
    
    ind_x = find(Hx_surf == 0);
    ind_y = find(Hy_surf == 0);
    ind_z = find(Hz_surf == 0);
    H_x_slice = HX(:,:,z_level);
    H_y_slice = HY(:,:,z_level);
    H_z_slice = HZ(:,:,z_level);
    Hx_surf(ind_x) = Hx_surf(ind_x) + H_x_slice(ind_x);
    Hy_surf(ind_y) = Hy_surf(ind_y) + H_y_slice(ind_y);
    Hz_surf(ind_z) = Hz_surf(ind_z) + H_z_slice(ind_z);
    
%     figure()
    quiver(X(:,:,1),Y(:,:,1),Hx_surf,Hy_surf);
    Hx_surf(141,171)
    Z(1,1,z_level)
    pause(0.5)
    
end

figure();
quiver(X(:,:,1),Y(:,:,1),Hx_surf,Hy_surf);


%% Calculate surface normal vector, take cross product with H field across surface to find surface current
% Surface normal vector = grad(f(x,y,z))/normalize

J_surf_tan_x = zeros(size(Hx_surf));
J_surf_tan_y = zeros(size(Hx_surf));

for i=1:length(X(:,1,1))
    for j=1:length(Y(1,:,1))
       normalVector = [-0.02*X(i,j,1),0,1]';
       unitVector = normalVector/norm(normalVector);
       Jsurf = cross(unitVector,[Hx_surf(i,j),Hy_surf(i,j),Hz_surf(i,j)]');
       J_surf_tan_x(i,j) = Jsurf(1)+Jsurf(3);
       J_surf_tan_y(i,j) = Jsurf(2);
    end
end

figure();
quiver(X(:,:,1),Y(:,:,1),Hx_surf,Hy_surf);
figure();
quiver(X(:,:,1),Y(:,:,1),J_surf_tan_x,J_surf_tan_y);


%% Calculate surface current density along x = -4.0 (~mid point of coil)

Xsurf = round(X(:,:,1),3,'significant');
Ysurf = round(Y(:,:,1),3,'significant');
%Find x = 0 index
zero_ind=find(Xsurf==-4.0);
% zero_ind = zero_ind(end);


%Make integral current go to zero:
JX_lineDensity = J_surf_tan_x(zero_ind)+abs(sum(J_surf_tan_x(zero_ind)))/length(J_surf_tan_x(zero_ind));
JY_lineDensity = J_surf_tan_y(zero_ind)+abs(sum(J_surf_tan_y(zero_ind)))/length(J_surf_tan_y(zero_ind));

figure();
plot(Ysurf(zero_ind),(JX_lineDensity));
hold on;
plot(Ysurf(zero_ind),(JY_lineDensity));

% figure();
% plot(Ysurf(zero_ind),cumsum(JX_lineDensity));
% 
% figure();
% plot(Ysurf(zero_ind),abs(cumsum(JX_lineDensity)));
% 



%% Calculate cross sections with equal currents

N = 4; %Number of compensation coil turns
I = 1; %Driving current
dx = (x(2)-x(1))*1e-3;
zci = @(v) find(diff(sign(v)));
ind=zci(JX_lineDensity);
ind1 = ind(1);
ind2=ind(2)+1;
zero_ind_const=176;

segment1 = fliplr(JX_lineDensity(1:ind1)'*dx); %Fliplr to start from the more current intense side of the curve (the side that matters more)
current1 = sum(segment1);
Itot1 = I*(N/2);
currentSegmentSum1 = Itot1/(N/2);
length1 = Xsurf(ind1,zero_ind_const)-Xsurf(1,zero_ind_const);

segment2 = JX_lineDensity(ind1+1:ind2-1)'*dx;
current2 = sum(segment2);
Itot2 = I*(N);
currentSegmentSum2 = Itot2/(N);
length2 = Xsurf(ind2-1,zero_ind_const)-Xsurf(ind1+1,zero_ind_const);


segment3 = JX_lineDensity(ind2:end)'*dx;
current3 = sum(segment3);
Itot3 = I*(N/2);
currentSegmentSum3 = Itot3/(N/2);
length3 = Xsurf(end,zero_ind_const)-Xsurf(ind2,zero_ind_const);

% Calculate optimal, even, number of compensation coil turns (assuming 1 A current):
optimalN = round_even((current2)/I)

subSegmentCurrents = [0,0,0];
i=1;
j=1;
k=1;
subStartInd = [i,j,k];
indices = {[],[],[]};
widths = {[],[],[]};
subSegmentAverages = {[],[],[]};

for m = 1:N
    if m <= N/2
    while subSegmentCurrents(1) < currentSegmentSum1
        subSegmentCurrents(1) = abs(sum(segment1(subStartInd(1):i)));
        i=i+1;
        if(i>length(segment1))
            i=i-1;
        break;
        end
    end
    end
    
    while subSegmentCurrents(2) < currentSegmentSum2
        subSegmentCurrents(2) = abs(sum(segment2(subStartInd(2):j)));
        j=j+1;
        if(j>length(segment2))
            j=j-1;
        break;
        end
    end
    
    if m > N/2
    while subSegmentCurrents(3) < currentSegmentSum3
        subSegmentCurrents(3) = abs(sum(segment1(subStartInd(3):k)));
        k=k+1;
        if(k>length(segment3))
            k=k-1;
        break;
        end
    end
    k=k-1;
    end
    

if m > N/2
    widths(2:3) = {[widths{2},Xsurf(j,zero_ind_const)-Xsurf(subStartInd(2),zero_ind_const)],[widths{3},Xsurf(k,zero_ind_const)-Xsurf(subStartInd(3),zero_ind_const)]};
    indices(2:3) = {[indices{2},j],[indices{3},k]};
    subSegmentAverages(2:3) = {[subSegmentAverages{2},sum(segment2(subStartInd(2):j))/length(segment2(subStartInd(2):j))], ...
    [subSegmentAverages{3},sum(segment3(subStartInd(3):k))/length(segment3(subStartInd(3):k))]};

else
    widths(1:2) = {[widths{1},Xsurf(i,zero_ind_const)-Xsurf(subStartInd(1),zero_ind_const)],[widths{2},Xsurf(j,zero_ind_const)-Xsurf(subStartInd(2),zero_ind_const)]};
    indices(1:2) = {[indices{1},i],[indices{2},j]};
    subSegmentAverages(1:2) = {[subSegmentAverages{1},sum(segment1(subStartInd(1):i))/length(segment1(subStartInd(1):i))], ...
    [subSegmentAverages{2},sum(segment2(subStartInd(2):j))/length(segment2(subStartInd(2):j))]};

end

subStartInd = [i+1,j+1,k+1];
subSegmentCurrents = [0,0,0];

    
end

% Flip back array
indices{1} = fliplr(length(segment1)-indices{1}+1);
segment1 = fliplr(segment1);

figure();
plot(Xsurf(:,zero_ind_const),(JX_lineDensity));
hold on;
scatter(Xsurf([indices{1}],zero_ind_const),zeros(size(Ysurf([indices{1}],zero_ind_const))));
scatter(Xsurf([indices{2}+ind1],zero_ind_const),zeros(size(Ysurf([indices{2}+ind1],zero_ind_const))));
scatter(Xsurf([indices{3}+ind2],zero_ind_const),zeros(size(Ysurf([indices{3}+ind2],zero_ind_const))));


%% Calculate streamline seed points and the corresponding streamlines

J_surf_abs = sqrt(J_surf_tan_x.^2 + J_surf_tan_y.^2);
J_Abs = J_surf_abs;

index_list = [ind1,indices{2}+ind1];
seedPointIndices=round(mean([index_list(1:end-1);index_list(2:end)]));


[max_value,idx]=max(J_Abs(:));
[xStart,yStart,z,w]=ind2sub(size(J_Abs),idx);
streamLineList = {};
streamLineLength = 15000;
streamLineDx = 2;

%Optional fix for streamlines that are too close to J = 0
% fixedY = Y;
% fixedY(seedPointIndices(1),zero_ind) = Y(seedPointIndices(1),zero_ind) + 0.5;
% fixedY(seedPointIndices(end),zero_ind) = Y(seedPointIndices(end),zero_ind) - 0.5;

for i = 1:length(seedPointIndices)
    startx = Xsurf(zero_ind_const,seedPointIndices(i));
    starty = Ysurf(zero_ind_const,seedPointIndices(i));
    XY = stream2(Xsurf,Ysurf,J_surf_tan_x,J_surf_tan_y,startx,starty,[streamLineDx,streamLineLength]);
    XY=cell2mat(XY);
    streamLineList{end+1} = XY;
    hold on;
end

%Crop streamlines to a single revolution:
for i = 1:length(streamLineList)
    streamLine = streamLineList{i};
    zci = @(v) find(diff(sign(v)));
    ind=zci(streamLine(:,1));
    streamLineList{i} = streamLine(1:ind(3)-10,:);
end

streamLineList

figure();
for i = 1:length(streamLineList)
    XY=streamLineList{i};
    plot(XY(:,1),XY(:,2));
    hold on;
end
