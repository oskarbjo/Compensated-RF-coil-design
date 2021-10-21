%% surface currents

% Read H fields at ground plane
close all;
load Hfield_data.mat;

x=Hfieldongroundplane41amp(:,1);
y=Hfieldongroundplane41amp(:,2);
Hx=Hfieldongroundplane41amp(:,3);
Hy=Hfieldongroundplane41amp(:,4);
Jx = [];
Jy = [];
zeroVector = zeros(length(Hx),1);
oneVector = ones(length(Hx),1);
J = cross([Hx,Hy,zeroVector],[zeroVector,zeroVector,oneVector]);
Jx = J(:,1);
Jy = J(:,2);
Jabs = sqrt(Jx.^2 + Jy.^2);


%calculate x y grid, assuming perfect square grid
[dummie,xmin_ind] = min(x);
[dummie,xmax_ind] = max(x);
[X,Y] = meshgrid(x(xmin_ind:xmax_ind),x(xmin_ind:xmax_ind));
JX = reshape(Jx,size(X))';
JY = reshape(Jy,size(X))';
J_Abs = sqrt(JX.^2 + JY.^2);


figure();
quiver(X,Y,JX,JY);


%% Calculate surface current density along x = 0


%Make integral current go to zero:
JX_lineDensity = JX(:,701)+abs(sum(JX(:,701)))/length(JX(:,701));
JY_lineDensity = JY(:,701)+abs(sum(JY(:,701)))/length(JY(:,701));

figure();
plot(Y(:,701),(JX_lineDensity));
hold on;
plot(Y(:,701),(JY_lineDensity));

figure();
plot(Y(:,701),cumsum(JX_lineDensity));

figure();
plot(Y(:,701),abs(cumsum(JX_lineDensity)));


%% Calculate cross sections with equal currents

N = 6; %Number of compensation coil turns
I = 1; %Driving current
dx = 0.05e-3;
zci = @(v) find(diff(sign(v)));
ind=zci(JX_lineDensity);
ind1 = ind(1);
ind2=ind(2)+1;

segment1 = fliplr(JX_lineDensity(1:ind1)'*dx); %Fliplr to start from the more current intense side of the curve (the side that matters more)
current1 = sum(segment1);
Itot1 = I*(N/2);
currentSegmentSum1 = Itot1/(N/2);
length1 = Y(ind1,701)-Y(1,701);

segment2 = JX_lineDensity(ind1+1:ind2-1)'*dx;
current2 = sum(segment2);
Itot2 = I*(N);
currentSegmentSum2 = -Itot2/(N);
length2 = Y(ind2-1,701)-Y(ind1+1,701);


segment3 = JX_lineDensity(ind2:end)'*dx;
current3 = sum(segment3);
Itot3 = I*(N/2);
currentSegmentSum3 = Itot3/(N/2);
length3 = Y(end,701)-Y(ind2,701);

%Calculate optimal, even, number of compensation coil turns (assuming 1 A current):
optimalN = round_even((-current2)/I)

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
        subSegmentCurrents(1) = sum(segment1(subStartInd(1):i));
        i=i+1;
        if(i>length(segment1))
            i=i-1;
        break;
        end
    end
    end
    
    while subSegmentCurrents(2) > currentSegmentSum2
        subSegmentCurrents(2) = sum(segment2(subStartInd(2):j));
        j=j+1;
        if(j>length(segment2))
            j=j-1;
        break;
        end
    end
    
    if m > N/2
    while subSegmentCurrents(3) < currentSegmentSum3
        subSegmentCurrents(3) = sum(segment1(subStartInd(3):k));
        k=k+1;
        if(k>length(segment3))
            k=k-1;
        break;
        end
    end
    k=k-1;
    end
    

if m > N/2
    widths(2:3) = {[widths{2},Y(j,701)-Y(subStartInd(2),701)],[widths{3},Y(k,701)-Y(subStartInd(3),701)]};
    indices(2:3) = {[indices{2},j],[indices{3},k]};
    subSegmentAverages(2:3) = {[subSegmentAverages{2},sum(segment2(subStartInd(2):j))/length(segment2(subStartInd(2):j))], ...
    [subSegmentAverages{3},sum(segment3(subStartInd(3):k))/length(segment3(subStartInd(3):k))]};

else
    widths(1:2) = {[widths{1},Y(i,701)-Y(subStartInd(1),701)],[widths{2},Y(j,701)-Y(subStartInd(2),701)]};
    indices(1:2) = {[indices{1},i],[indices{2},j]};
    subSegmentAverages(1:2) = {[subSegmentAverages{1},sum(segment1(subStartInd(1):i))/length(segment1(subStartInd(1):i))], ...
    [subSegmentAverages{2},sum(segment2(subStartInd(2):j))/length(segment2(subStartInd(2):j))]};

end

subStartInd = [i+1,j+1,k+1];
subSegmentCurrents = [0,0,0];

    
end

%Flip back array
indices{1} = fliplr(length(segment1)-indices{1}+1);
segment1 = fliplr(segment1);

figure();
plot(Y(:,701),(JX_lineDensity));
hold on;
scatter(Y([indices{1}],701),zeros(size(Y([indices{1}],701))));
scatter(Y([indices{2}+ind1],701),zeros(size(Y([indices{2}+ind1],701))));
scatter(Y([indices{3}+ind2],701),zeros(size(Y([indices{3}+ind2],701))));



%% Calculate streamline seed points and the corresponding streamlines


index_list = [ind1,indices{2}+ind1];
seedPointIndices=round(mean([index_list(1:end-1);index_list(2:end)]));


[max_value,idx]=max(J_Abs(:));
[xStart,yStart,z,w]=ind2sub(size(J_Abs),idx);
streamLineList = {};
streamLineLength = 15000;
streamLineDx = 2;

%Optional fix for streamlines that are too close to J = 0
fixedY = Y;
fixedY(seedPointIndices(1),701) = Y(seedPointIndices(1),701) + 0.5;
fixedY(seedPointIndices(end),701) = Y(seedPointIndices(end),701) - 0.5;

for i = 1:length(seedPointIndices)
    startx = X(seedPointIndices(i),701);
    starty = fixedY(seedPointIndices(i),701);
    XY = stream2(X,Y,JX,JY,startx,starty,[streamLineDx,streamLineLength]);
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

maxWidth = 8; %[mm]
Npts = 10001; %Should be an odd number
pts = linspace(-maxWidth/2,maxWidth/2,Npts);
dl = pts(2)-pts(1);
traceOutlines = cell(length(streamLineList),1);

for i = 1:length(traceOutlines)
    traceOutlines{i} = zeros(length(streamLineList{i}),4);
end
figure();
for i = 1:length(streamLineList)
    streamLine = streamLineList{i};
    for j = 1:length(streamLine)
        j
        % 1. Calculate line perpendicular to streamline
        x=streamLine(j,1);
        y=streamLine(j,2);
        Jx_interp = interp2(X,Y,JX,x,y); %Jx and Jy on streamline
        Jy_interp = interp2(X,Y,JY,x,y);
        unitVector = [Jx_interp,Jy_interp]/sqrt(Jx_interp^2+Jy_interp^2);
        perpVector = [unitVector(2),-unitVector(1)];
        integrationPathX = x + pts*perpVector(1);
        integrationPathY = y + pts*perpVector(2);
        if j == 350
            j
            plot(streamLine(:,1),streamLine(:,2));
            hold on;
            plot(integrationPathX,integrationPathY);
        end
        
        % 2. Interpolate Jx and Jy along the perpendicular line
        Jx_line = interp2(X,Y,JX,integrationPathX,integrationPathY);
        Jy_line = interp2(X,Y,JY,integrationPathX,integrationPathY);
        
        % 3. Take the dot product between J_line and unitVector 
        a=[Jx_line',Jy_line'];
        b=unitVector';
        dotProduct = a*b;
        dotProduct = dotProduct';
        % 4. and 5. Perform line integral. First starting from zero and integrating 
        % in the one direction, then starting from zero and integrating in 
        % the other direction until each integral = I/2
        integralDir1 = 0;
        integralDir2 = 0;
        startingPoint = (numel(dotProduct)+1)/2;
        m = 1;
        n = 1;
        flippedDotProduct = fliplr(dotProduct);
        while integralDir1 < (I/N)
            integralDir1 = sum(flippedDotProduct(startingPoint:startingPoint+m))*dl*0.001;
            m=m+1;
            if(startingPoint+m > length(dotProduct))
                m=m-1;
                break;
            elseif integralDir1 < 0
                    break;
            elseif(i==1 || i==length(streamLine)) %optional if statement that fixes the smallest traces
                if integralDir1 > 0.8* I/N
                    break;
                end
            end
        end
        while integralDir2 < (I/N)
            integralDir2 = sum(dotProduct(startingPoint:startingPoint+n))*dl*0.001;
            n=n+1;
            if(startingPoint+n > length(dotProduct))
                n=n-1;
                break;
            elseif integralDir2 < 0
                    break;
            elseif(i==1 || i==length(streamLine)) %optional if statement that fixes the smallest traces
                if integralDir2 > 0.8* I/N
                    break;
                end
            end
        end
        %6. Save coordinates of trace outline
        traceOutlines{i}(j,:) = [integrationPathX(startingPoint-m),integrationPathY(startingPoint-m),integrationPathX(startingPoint+n),integrationPathY(startingPoint+n)];
        
    end
    
    color=[rand,rand,rand];
    plot(traceOutlines{i}(:,1),traceOutlines{i}(:,2),'color',color,'linewidth',2);
    hold on;
    plot(traceOutlines{i}(:,3),traceOutlines{i}(:,4),'color',color,'linewidth',2);
    XY=streamLineList{i};
    plot(XY(:,1),XY(:,2),'black','linewidth',2);
    drawnow;
end

quiver(X,Y,JX,JY);

%% Concatenate data to be exported as a single polygon

output=[];
finalpoint=[];
for i = 1:length(traceOutlines)
   
    output = [output;traceOutlines{i}(:,1:2)];
    
end
for i = length(traceOutlines):-1:1

    output = [output;fliplr(traceOutlines{i}(:,3:4)')'];
    
end
output = [output;finalpoint];
output = [output(1,:);output];