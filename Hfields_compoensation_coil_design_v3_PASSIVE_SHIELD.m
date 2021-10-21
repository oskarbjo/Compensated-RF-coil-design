%% Import H field data and calculate surface currents

clc;
% close all;
% Use 3D volume field export for Z=Ground plane
% path="C:\Users\objoerkqvist\OneDrive - ETH Zurich\Documents\CST\Harvesting coil designs\vertical shielded coil v4\GP_currents_v4.txt";
path="C:\Users\objoerkqvist\OneDrive - ETH Zurich\Documents\CST\Harvesting coil designs\vertical shielded coil v4\shield_currents_check_6e-11_v4.txt";
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
Jx = J(:,1);
Jy = J(:,2);

%calculate x y grid, assuming perfect square grid
[dummie,xmin_ind] = min(x);
[dummie,xmax_ind] = max(x);
shift = max(x)-max(y);
[X,Y] = meshgrid(x(xmin_ind:xmax_ind),x(xmin_ind:xmax_ind)-shift);
JX = reshape(Jx,size(X))';
JY = reshape(Jy,size(X))';
J_Abs = sqrt(JX.^2 + JY.^2);

figure();
h1=quiver(X,Y,JX,JY);
set(h1,'AutoScale','on', 'AutoScaleFactor', 2);

%% Calculate surface current density along x = 0


%Find x = 0 index
zero_ind=find(X==0)/length(X);
zero_ind = zero_ind(end); 

%Make integral current go to zero:
JX_lineDensity = JX(:,zero_ind)+abs(sum(JX(:,zero_ind)))/length(JX(:,zero_ind));
JY_lineDensity = JY(:,zero_ind)+abs(sum(JY(:,zero_ind)))/length(JY(:,zero_ind));

figure();
plot(Y(:,zero_ind),(JX_lineDensity));
xlabel('x position [mm]')
ylabel('Surface current [A/m]')
hold on;
% plot(Y(:,zero_ind),(JY_lineDensity));

figure();
plot(Y(:,zero_ind),cumsum(JX_lineDensity));

figure();
plot(Y(:,zero_ind),abs(cumsum(JX_lineDensity)));


%% Divide middle section into N segments of equal current

N = 16; %Number of divisions (should be even number)
dx = (x(2)-x(1))*1e-3;
zci = @(v) find(diff(sign(v))); %zero crossing index finder function
ind=zci(JX_lineDensity);
minDistanceFromZC = 0; %*** OPTIONAL *** Move shields away from local zero 
ind1 = ind(1) + minDistanceFromZC;
ind2=ind(2)+1 - minDistanceFromZC;


segment2 = JX_lineDensity(ind1+1:ind2-1)'*dx;
current2 = sum(segment2);
length2 = Y(ind2-1,zero_ind)-Y(ind1+1,zero_ind);


subSegmentCurrents = 0;
j=1;

subStartInd = j;
indices = [];
widths = [];
subSegmentAverages = [];
seedPointInd=[];
findSeedPoint=1;
for m = 1:N
    while subSegmentCurrents < abs(current2/N)
        subSegmentCurrents = abs(sum(segment2(subStartInd:j)));
        j=j+1;
        if(subSegmentCurrents > abs(current2)/N/2 && findSeedPoint==1)
            seedPointInd=[seedPointInd,j];
            findSeedPoint=0;
        end
        if(j>length(segment2))
            j=j-1;
        break;
        end
    end
    
findSeedPoint=1;    
indices(m) = j;
subStartInd = j+1;
subSegmentCurrents = 0;

    
end


figure();
plot(Y(:,zero_ind),(JX_lineDensity));
hold on;
scatter(Y([indices+ind1],zero_ind),zeros(size(Y([indices+ind1],zero_ind))));

[a1,b1]=unique(Y([indices+ind1]));
b1(end)
indices=indices(1:b1(end));
I=2/numel(indices);
%% Calculate streamline seed points and the corresponding streamlines
%%%%%%%%%%%%%%%%%%%%%
% Option 1: Midway between calculated index points
index_list = [ind1,indices+ind1];
seedPointIndices=round(mean([index_list(1:end-1);index_list(2:end)]));
shiftPts = [-1,0,-3,-2,-1,1,0,0,-4,-5,0,0,0,0];
seedPointIndices = seedPointIndices - shiftPts;
% seedPointInd(:)=0;
%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%
% Option 2: Point that splits current in half
% shiftPts1=0;
% shiftPts2=0;
% shiftPts3=0;
% shiftPts4=0;
% shiftPts5=0;
% shiftPts6=0;
% shiftPts7=0;
% shiftPts8=0;
% shiftPts9=0;
% shiftPts10=0;
% shiftPts = zeros(1,numel(seedPointInd));
% seedPointIndices = ind1 + seedPointInd - shiftPts;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



[max_value,idx]=max(J_Abs(:));
[xStart,yStart,z,w]=ind2sub(size(J_Abs),idx);
streamLineList = {};
streamLineLength = 15000;
streamLineDx = [6,5,5,5,4,4,3]; %Step size for each individual shield trace - should have N/2 elements, where N is number of shield conductors
streamLineDx = [fliplr(streamLineDx),streamLineDx];
cutPoints = [2,2,2,2,2,2,4];%Cutpoints * streamLineDx = length cut from streamline to facilitate interconnections
cutPoints = [fliplr(cutPoints),cutPoints];

%Optional fix for streamlines that are too close to J = 0
fixedY = Y;
% fixedY(seedPointIndices(1),zero_ind) = Y(seedPointIndices(1),zero_ind) + 0.4;
% fixedY(seedPointIndices(end),zero_ind) = Y(seedPointIndices(end),zero_ind) - 0.4;

xOffsets = linspace(0,6,numel(seedPointIndices)/2); %For shifting seedpoints forward to facilitate interconnections
xOffsets = [xOffsets, fliplr(xOffsets)];

for i = 1:length(seedPointIndices)
    startx = X(seedPointIndices(i),zero_ind)-xOffsets(i);
    starty = fixedY(seedPointIndices(i),zero_ind);
    XY = stream2(X,Y,JX,JY,startx,starty,[streamLineDx(i),streamLineLength]);
    XY=cell2mat(XY);
    streamLineList{end+1} = XY;
end


%Crop streamlines to a single revolution:
for i = 1:length(streamLineList)
    streamLine = streamLineList{i};
    zci = @(v) find(diff(sign(v)));
    ind=zci(streamLine(:,1)+xOffsets(i));
    streamLineList{i} = streamLine(1:ind(3)-cutPoints(i),:);
end

%%%%%%% Choose streamlines such that their pairs have the same size
targetAreas = [];
for j = 1:length(seedPointIndices)/2
   targetAreas(j) = polyarea(streamLineList{j}(:,1),streamLineList{j}(:,2));
    
end
targetAreas = fliplr(targetAreas);

correctLengthsBoolean = ones(length(seedPointIndices)/2,1);
shifts=zeros(length(seedPointIndices)/2,1);

while sum(correctLengthsBoolean ~= 0)
    
for i = (length(streamLineList)/2+1):length(streamLineList)
    jj=i-length(streamLineList)/2;
    shieldArea = polyarea(streamLineList{i}(:,1),streamLineList{i}(:,2))
    targetAreas(jj)
    shiftSign=1;
    if shieldArea < targetAreas(jj)
        shifts(jj) = shifts(jj)+1;
        startx = X(seedPointIndices(i),zero_ind)-xOffsets(i);
        starty = fixedY(seedPointIndices(i)-shiftSign*shifts(jj),zero_ind);
        XY = stream2(X,Y,JX,JY,startx,starty,[streamLineDx(i),streamLineLength]);
        XY=cell2mat(XY);
        streamLineList{i} = XY;
        streamLine = streamLineList{i};
        zci = @(v) find(diff(sign(v)));
        ind=zci(streamLine(:,1)+xOffsets(i));
        pastStreamLines = streamLineList;
        streamLineList{i} = streamLine(1:ind(3)-cutPoints(i),:);
        if abs(1-polyarea(pastStreamLines{i}(:,1),pastStreamLines{i}(:,2))/targetAreas(jj)) > abs(1-polyarea(pastStreamLines{i}(:,1),pastStreamLines{i}(:,2))/targetAreas(jj))
            streamLineList{i}=pastStreamLines{i};
        end

    else
        correctLengthsBoolean(jj)=0;
    end
    
    
end
allAreas = [];
for j = 1:length(seedPointIndices)
   allAreas(j) = polyarea(streamLineList{j}(:,1),streamLineList{j}(:,2)); 
end


end

streamLineList

figure();
for i = 1:length(streamLineList)
    XY=streamLineList{i};
    plot(XY(:,1),XY(:,2),'linewidth',2);
    hold on;
end


streamLineLengths = [];
areaList=[];

for i = 1:length(streamLineList)
    len = 0;
    for j = 1:length(streamLineList{i})-1
        len = len + sqrt( (streamLineList{i}(j+1,1)-streamLineList{i}(j,1))^2 + (streamLineList{i}(j+1,2)-streamLineList{i}(j,2))^2 );
    end
    areaList=[areaList;polyarea(streamLineList{i}(:,1),streamLineList{i}(:,2))];
    streamLineLengths = [streamLineLengths,len];
end

streamLineLengths
allAreas./(fliplr(allAreas))



%% Calculate surface current integrals

% Process:
% 1. Calculate a line perpendicular to streamline to integrate surface 
% currents along (p/m L/2). The length of this line is the max trace width.
% 2. Interpolate Jx and Jy components along this line.
% 3. Take the dot product between [Jx,Jy] and a unit vector that is 
% parallel with the streamline
% 4. Integrate this dot product from its mid point in the one direction
% until the integral = I/2
% 5. Integrate this dot product from its mid point in the other direction
% until the integral = I/2
% 6. Lower and upper bounds of integral give width of trace at this point
% on the streamline

maxWidth = 0.15; %[mm], adjust widths according to each trace length
traceStandardWidth=maxWidth(1);
Npts = 3; %Should be an odd number
traceOutlines = cell(length(streamLineList),1);
splitLines={{[]}};
splitWidthLimit=0.5; %Limit for splitting lines
currentLimit=I; %I/N


for i = 1:length(traceOutlines)
    traceOutlines{i} = zeros(length(streamLineList{i}),4);
end
figure();
for i = 1:length(streamLineList)
    streamLine = streamLineList{i};
    % Part of step 1: Calculate all interpolated surface current values at once to
    % make code run faster
    xVector=streamLine(:,1);
    yVector=streamLine(:,2);
    Jx_interp_vector = interp2(X,Y,JX,xVector,yVector);
    Jy_interp_vector = interp2(X,Y,JY,xVector,yVector);
    pts = linspace(-maxWidth/2,maxWidth/2,Npts);
    dl=pts(2)-pts(1);
    for j = 1:length(streamLine)
        j
        '1.'
        tic
        % 1. Calculate line perpendicular to streamline
        unitVector = [Jx_interp_vector(j),Jy_interp_vector(j)]/sqrt(Jx_interp_vector(j)^2+Jy_interp_vector(j)^2);
        perpVector = [unitVector(2),-unitVector(1)];
        integrationPathX = xVector(j) + pts*perpVector(1);
        integrationPathY = yVector(j) + pts*perpVector(2);
        if mod(j,15)==0
            if j==2
            plot(streamLine(:,1),streamLine(:,2));
            end
            hold on;
            plot(integrationPathX,integrationPathY);
            axis equal;
            drawnow;
        end
        x1=xVector(j)+perpVector(1)*traceStandardWidth/2;
        x2=xVector(j)-perpVector(1)*traceStandardWidth/2;
        y1=yVector(j)+perpVector(2)*traceStandardWidth/2;
        y2=yVector(j)-perpVector(2)*traceStandardWidth/2;
        traceOutlines{i}(j,:) = [x1,y1,x2,y2];
        toc
        
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
    [m,n]=max(abs(streamLineList{i}(:,2)));
    o1=0;
    while abs(streamLineList{i}(n+o1,1)-streamLineList{i}(n,1)) < capLength/2
        o1=o1+1;
    end
    o2=0;
    while abs(streamLineList{i}(n-o2,1)-streamLineList{i}(n,1)) < capLength/2
        o2=o2+1;
    end
    capacitorTraceoutlineList{i*2-1}=traceOutlines{i}(1:n-o1,:);
    capacitorTraceoutlineList{i*2}=traceOutlines{i}(n+o2:end,:);
    
    
end

figure();
for i=1:length(capacitorTraceoutlineList)
    plot(capacitorTraceoutlineList{i}(:,1),capacitorTraceoutlineList{i}(:,2),'color','b','linewidth',1);
    hold on;
    plot(capacitorTraceoutlineList{i}(:,3),capacitorTraceoutlineList{i}(:,4),'color','b','linewidth',1);
end


%% Find points to place wires and capacitors on

wireMountPoints = [];
capacitorMountPoints = [];

for i = 1:length(traceOutlines)
    wireMountPoints(i,1:2) = [capacitorTraceoutlineList{i*2-1}(1,3)-(capacitorTraceoutlineList{i*2-1}(1,3)-capacitorTraceoutlineList{i*2-1}(1,1))/2, capacitorTraceoutlineList{i*2-1}(1,4)-(capacitorTraceoutlineList{i*2-1}(1,4)-capacitorTraceoutlineList{i*2-1}(1,2))/2];
    wireMountPoints(length(traceOutlines)-i+1,3:4) = [capacitorTraceoutlineList{i*2}(end,3)-(capacitorTraceoutlineList{i*2}(end,3)-capacitorTraceoutlineList{i*2}(end,1))/2, capacitorTraceoutlineList{i*2}(end,4)-(capacitorTraceoutlineList{i*2}(end,4)-capacitorTraceoutlineList{i*2}(end,2))/2];
    capacitorMountPoints(i,1:2) = [capacitorTraceoutlineList{i*2-1}(end,3)-(capacitorTraceoutlineList{i*2-1}(end,3)-capacitorTraceoutlineList{i*2-1}(end,1))/2, capacitorTraceoutlineList{i*2-1}(end,4)-(capacitorTraceoutlineList{i*2-1}(end,4)-capacitorTraceoutlineList{i*2-1}(end,2))/2];
    capacitorMountPoints(i,3:4) = [capacitorTraceoutlineList{i*2}(1,3)-(capacitorTraceoutlineList{i*2}(1,3)-capacitorTraceoutlineList{i*2}(1,1))/2, capacitorTraceoutlineList{i*2}(1,4)-(capacitorTraceoutlineList{i*2}(1,4)-capacitorTraceoutlineList{i*2}(1,2))/2];
end
scatter(wireMountPoints(:,1),wireMountPoints(:,2),'r');
scatter(wireMountPoints(:,3),wireMountPoints(:,4),'r');
scatter(capacitorMountPoints(:,1),capacitorMountPoints(:,2),'b');
scatter(capacitorMountPoints(:,3),capacitorMountPoints(:,4),'g');

%% Write polygons and conncetion points to file

path = "C:\Users\objoerkqvist\OneDrive - ETH Zurich\Documents\MATLAB\Coil design\polygon_export\";

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

path="C:\Users\objoerkqvist\OneDrive - ETH Zurich\Documents\MATLAB\Coil design\wire_cap_points\";
fileID = fopen(path + "wirePoints" + '.txt','w');
fprintf(fileID,'%d %d %d %d \n',wireMountPoints');
fclose(fileID);
fileID = fopen(path + "capPoints" + '.txt','w');
fprintf(fileID,'%d %d %d %d \n',capacitorMountPoints');
fclose(fileID);
