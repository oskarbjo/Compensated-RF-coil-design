% SHIELD DESIGN - NEW SIMPLIFIED VERSION



%% Import H field data and calculate surface currents

%Remember that it's easier to mirror one side of the shield to the other,
%rather than trying to make the lobes on each side the exact same size
%manually. This works just as well!

clc;
close all;
clear all;
% Use 3D volume field export for Z=Ground plane
% path="C:\Users\objoerkqvist\OneDrive - ETH Zurich\Documents\CST\Harvesting coil designs\Vertical shielded coil v7 quadrature\POS1_Hfield.txt"; %Current shield!
path="C:\Users\objoerkqvist\OneDrive - ETH Zurich\Documents\CST\Harvesting coil designs\Harvest_coil_v_2_0\H field export\closest_coil.txt";

data = importdata(path);
data=struct2cell(data);


x_ind = 1;
y_ind = 2;
z_ind = 3;
Hx_ind = 4;
Hy_ind = 6;
Hz_ind = 8;

% XY should span ground plane
x = data{1}(:,x_ind);
y = data{1}(:,y_ind);
z = data{1}(:,z_ind);
Hx = data{1}(:,Hx_ind);
Hy = data{1}(:,Hy_ind);
Hz = data{1}(:,Hz_ind);


zeroVector = zeros(length(Hx),1);
oneVector = ones(length(Hx),1);
J = cross([Hx,Hy,Hz],[zeroVector,zeroVector,oneVector]);
Jx = J(:,2);
Jy = J(:,1);



%calculate x y grid, assuming perfect square grid

[dummie,xmin_ind] = min(x);
[dummie,xmax_ind] = max(x);
shift = max(x)-max(y);
[X,Y] = meshgrid(x(xmin_ind:xmax_ind),x(xmin_ind:xmax_ind)-shift);
Y=round(Y,5);
X=round(X,5);

JX = reshape(Jx,size(X));
JY = reshape(Jy,size(X));

%Coil windings should be parallel to y - use following if needed
JXtemp = JX;
JXtmp = rot90(JX);
JY = rot90(JY);
JX=-JY;
JY=JXtmp;


J_Abs = sqrt(JX.^2 + JY.^2);


%ENTER SHIFTS (specific for current dataset) - set origin in middle of
%coils
% X=X-5; 
X = X - 7.5;
    

figure();
h1=quiver(X,Y,JX,JY);
set(h1,'AutoScale','on', 'AutoScaleFactor', 2);


%% Calculate surface current density along x = 0


%Coil conductors should be parallel to y. Seed points will be chosen in
%outer region (not directly under conductors)

%ENTER EXTENT OF SHIELD PCB
yextent = 75;
xextent = 75;

%ENTER NUMBER OF SHIELD CONDUCTORS (for one side, total will be x2):
N_cond = 15;

%ENTER NUMBER OF EDGE SAMPLES THAT SHOULD BE IGNORED:
edge = 2;


%Find 0 index
[rowX,colX,dat]=find(Y==0); %Search along y
zero_ind = rowX(end); 

[rowX,colX,dat]=find(X>xextent/2); %Search along x 
extent_ind = colX(1);

zci = @(v) find(diff(sign(v))); %zero crossing index finder function
zeroCrInd = zci(JY(zero_ind,:));
zeroCrInd = zeroCrInd(end);

current_density_cross_section = JY(zero_ind,zeroCrInd+edge:extent_ind-edge);
current_density_cross_section_pos = Y(zero_ind,zeroCrInd+edge:extent_ind-edge);

% Split the current cross section into N_cond equal parts:
current_share = abs(sum(current_density_cross_section)/N_cond);
index1=1;
index2=1;
seedPointIndices = zeros(N_cond,1);
for i = 1:N_cond
    try
        while(abs(sum(current_density_cross_section(index1:index2))) < current_share)
            index2=index2+1;
        end
    catch
    end
    index1=index2;
    seedPointIndices(i) = edge+zeroCrInd+round(mean([index1,index2]));
end


figure()
plot(Y(:,1),JY(zero_ind,:))
hold on
scatter(Y(seedPointIndices,1),zeros(numel(X(seedPointIndices,zero_ind)),1));




%% CALCULATE STREAMLINES

%ENTER STREAMLINE LENGTH
streamLineLength = 2000;


%ENTER STREAMLINE DX
streamLineDx=1.5;




streamLineList = {};
figure();
%Calculate raw streamlines
for i = 1:length(seedPointIndices)
    startx = X(zero_ind,seedPointIndices(i));
    starty = Y(zero_ind,seedPointIndices(i));
    XY = stream2(X,Y,JX,JY,startx,starty,[streamLineDx,streamLineLength]);
    XY=cell2mat(XY);
    streamLineList{end+1} = XY;
end

%Crop streamlines to one revolution:
for i = 1:length(streamLineList)
    streamLine = streamLineList{i};
    zci = @(v) find(diff(sign(v)));
    ind=zci(streamLine(:,2));
    streamLineList{i} = [[flipud(streamLine(1:ind(2),1)); streamLine(2:ind(2),1)], [flipud(streamLine(1:ind(2),2)); -streamLine(2:ind(2),2)]]; %Save cropped streamline and put center point first
end

%Get rid of points that overlap into negative x if applicable
for i = 1:length(streamLineList)
    streamLine = streamLineList{i};
    for j = 1:length(streamLine)
        if(streamLine(j,1) < 0.075)
            streamLine(j,1) = 0.075;
        end
    end
    streamLineList{i} = streamLine;
    
end

%Cut starting points for connectivity
cutDist = 1.5; %mm
offset=2; %EDIT!
step = 2; %EDIT!
startInd = offset:step:step*numel(streamLineList)+offset;
for i = 1:length(streamLineList)
    ind=startInd(i)+round(50/streamLineDx);
    dist=1e9;
    while(dist > cutDist)
        ind = ind+1;
        if(ind > numel(streamLineList{i}(:,1)))
            ind = 1;
        end
        dist = norm(streamLineList{i}(startInd(i),:) - streamLineList{i}(ind,:));        
    end
    if(ind <= numel(streamLineList{i}(:,1)) && ind > startInd(i))
        streamLineList{i} = [streamLineList{i}(startInd(i):ind,:)];
    elseif(ind < startInd(i))
        streamLineList{i} = [streamLineList{i}(startInd(i):end,:); streamLineList{i}(1:ind,:)];
    end
    plot(streamLineList{i}(:,1),streamLineList{i}(:,2));
    hold on;
    drawnow;
end


%% Move streamlines if applicable (might not be the case!):

movex = fliplr([-1.6, -1.3,-1.2,-1,-0.8,-0.7,-0.5,-0.4,-0.3,-0.2,-0.1,0,0,0,0]);
figure;
for i = 1:length(streamLineList)
    streamLineList{i}(:,1) = streamLineList{i}(:,1) + movex(i); 
    plot(streamLineList{i}(:,1),streamLineList{i}(:,2));
    hold on;
    drawnow;
end

%% Create outlines for conductors


width = 0.1; %[mm]
Npts = 3; %Should be an odd number
traceOutlines = cell(length(streamLineList),1);
figure;


traceOutlines={};

figure();
for i = 1:length(streamLineList)
    streamLine = streamLineList{i};
    xVector=streamLine(:,1);
    yVector=streamLine(:,2);
    for j = 1:length(streamLine)-1
        % 1. Calculate line perpendicular to streamline
        dx=xVector(j+1)-xVector(j);
        dy=yVector(j+1)-yVector(j);
        unitVector = [dx,dy]/sqrt(dx^2+dy^2);
        perpVector = [unitVector(2),-unitVector(1)];

        x1=xVector(j)+perpVector(1)*width/2;
        x2=xVector(j)-perpVector(1)*width/2;
        y1=yVector(j)+perpVector(2)*width/2;
        y2=yVector(j)-perpVector(2)*width/2;
        traceOutlines{i}(j,:) = [x1,y1,x2,y2];
    end

    color=[rand,rand,rand];
    plot(traceOutlines{i}(:,1),traceOutlines{i}(:,2),'color','b','linewidth',1);
    hold on;
    plot(traceOutlines{i}(:,3),traceOutlines{i}(:,4),'color','b','linewidth',1);
    XY=streamLineList{i};
    plot(XY(:,1),XY(:,2),'black','linewidth',1);
    drawnow;

end



%% Function that chops up the shields to make room for capacitors
capacitorTraceoutlineList = cell(2*length(streamLineList),1);
capLength = 1.2;
for i = 1:length(streamLineList)
    [m,n]=max(abs(streamLineList{i}(:,1)));
    o1=0;
    while abs(streamLineList{i}(n+o1,2)-streamLineList{i}(n,2)) < capLength/2
        o1=o1+1;
    end
    o2=0;
    while abs(streamLineList{i}(n-o2,2)-streamLineList{i}(n,2)) < capLength/2
        o2=o2+1;
    end

    capacitorTraceoutlineList{i*2-1}=traceOutlines{i}(1:n-o1,:);
    capacitorTraceoutlineList{i*2}=traceOutlines{i}(n+o2:end,:);
    
end

N=length(capacitorTraceoutlineList);
for i = 1:2:N
    capacitorTraceoutlineList{end+1} = [-capacitorTraceoutlineList{i}(:,1), capacitorTraceoutlineList{i}(:,2), -capacitorTraceoutlineList{i}(:,3), capacitorTraceoutlineList{i}(:,4)];
    capacitorTraceoutlineList{end+1} = [-capacitorTraceoutlineList{i+1}(:,1), capacitorTraceoutlineList{i+1}(:,2), -capacitorTraceoutlineList{i+1}(:,3), capacitorTraceoutlineList{i+1}(:,4)];
end

figure();
for i=1:length(capacitorTraceoutlineList)
    plot(capacitorTraceoutlineList{i}(:,1),capacitorTraceoutlineList{i}(:,2),'color','b','linewidth',1);
    hold on;
    plot(capacitorTraceoutlineList{i}(:,3),capacitorTraceoutlineList{i}(:,4),'color','b','linewidth',1);
end

%% Find points to place wires and capacitors on

N=length(capacitorTraceoutlineList)/2;
wireMountPoints = zeros(N,4);
capacitorMountPoints = zeros(N,4);
j=1;
for i = 1:N/2 %this is implemented in an absolutely stupid way
    x1=(capacitorTraceoutlineList{i*2-1}(1,1)+capacitorTraceoutlineList{i*2-1}(1,3))/2;
    x2=(capacitorTraceoutlineList{N+i*2}(end,1)+capacitorTraceoutlineList{N+i*2}(end,3))/2;
    y1=(capacitorTraceoutlineList{i*2-1}(1,2)+capacitorTraceoutlineList{i*2-1}(1,4))/2;
    y2=(capacitorTraceoutlineList{N+i*2}(end,2)+capacitorTraceoutlineList{N+i*2}(end,4))/2;
    wireMountPoints(j,:) = [x1,y1,x2,y2];
    j=j+1;
    x1=(capacitorTraceoutlineList{i*2}(end,1)+capacitorTraceoutlineList{i*2}(end,3))/2;
    x2=(capacitorTraceoutlineList{N+i*2-1}(1,1)+capacitorTraceoutlineList{N+i*2-1}(1,3))/2;
    y1=(capacitorTraceoutlineList{i*2}(end,2)+capacitorTraceoutlineList{i*2}(end,4))/2;
    y2=(capacitorTraceoutlineList{N+i*2-1}(1,2)+capacitorTraceoutlineList{N+i*2-1}(1,4))/2;
    wireMountPoints(j,:) = [x1,y1,x2,y2];
    j=j+1;
end

for i=1:N

    x1=(capacitorTraceoutlineList{i*2-1}(end,1)+capacitorTraceoutlineList{i*2-1}(end,3))/2;
    x2=(capacitorTraceoutlineList{i*2}(1,1)+capacitorTraceoutlineList{i*2}(1,3))/2;
    y1=(capacitorTraceoutlineList{i*2-1}(end,2)+capacitorTraceoutlineList{i*2-1}(end,4))/2;
    y2=(capacitorTraceoutlineList{i*2}(1,2)+capacitorTraceoutlineList{i*2}(1,4))/2;
    
    capacitorMountPoints(i,:) = [x1,y1,x2,y2];

end

scatter(wireMountPoints(:,1),wireMountPoints(:,2),'r');
hold on;
scatter(wireMountPoints(:,3),wireMountPoints(:,4),'r');
scatter(capacitorMountPoints(:,1),capacitorMountPoints(:,2),'b');
scatter(capacitorMountPoints(:,3),capacitorMountPoints(:,4),'g');

%% Write polygons and connection points to file

%Save the old files before new exports or the older CST files will not be
%generated correctly. Also delete the files in the write directory to
%ensure that no old polygons remain

path = "C:\Users\objoerkqvist\OneDrive - ETH Zurich\Documents\CST\Harvesting coil designs\Harvest_coil_v_2_0\shield data\polygons\";

for i = 1:length(capacitorTraceoutlineList)
    output = [];
    output = [output;capacitorTraceoutlineList{i}(:,1:2)];
    output = [output;fliplr(capacitorTraceoutlineList{i}(:,3:4)')'];
    output = [output;output(1,:)];
%     dlmwrite(path + "polygon" + num2str(i) + '.txt',output,'delimiter',' ');
    fileID = fopen(path + "polygon" + num2str(i) + '.txt','w');
    fprintf(fileID,'%d %d \n',output');
    fclose(fileID);
end

path="C:\Users\objoerkqvist\OneDrive - ETH Zurich\Documents\CST\Harvesting coil designs\Harvest_coil_v_2_0\shield data\wire_cap_points\";
fileID = fopen(path + "wirePoints" + '.txt','w');
fprintf(fileID,'%d %d %d %d \n',wireMountPoints');
fclose(fileID);
fileID = fopen(path + "capPoints" + '.txt','w');
fprintf(fileID,'%d %d %d %d \n',capacitorMountPoints');
fclose(fileID);