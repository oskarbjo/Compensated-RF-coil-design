%% Import H field data and calculate surface currents

clc;
close all;

path="C:\Users\objoerkqvist\OneDrive - ETH Zurich\Documents\CST\Harvesting coil designs\vertical shielded coil v3\surface_H_field_GPgap1p2_Nturns4_blend1_H10_L25_windinggap2p1.txt";
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
Jz = J(:,3);
Jabs = sqrt(Jx.^2 + Jy.^2 + Jz.^2);

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
set(h1,'AutoScale','on', 'AutoScaleFactor', 2)

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


%% Calculate cross sections with equal currents

N = 21; %Number of compensation coil turns (left + right side, should be even number)
I = 2/N; %Driving current (this equals the current flowing in each compensation winding)
dx = (x(2)-x(1))*1e-3;
zci = @(v) find(diff(sign(v))); %zero crossing index finder function
ind=zci(JX_lineDensity);
ind1 = ind(1);
ind2=ind(2)+1;

segment1 = fliplr(JX_lineDensity(1:ind1)'*dx); %Fliplr to start from the more current intense side of the curve (the side that matters more)
current1 = sum(segment1);
Itot1 = round(current1); %Current that SHOULD flow on the one side of the main coil
currentSegmentSum1 = Itot1/(N/2);
length1 = Y(ind1,zero_ind)-Y(1,zero_ind);

segment2 = JX_lineDensity(ind1+1:ind2-1)'*dx;
current2 = sum(segment2);
Itot2 = round(current2);
currentSegmentSum2 = Itot2/(N);
length2 = Y(ind2-1,zero_ind)-Y(ind1+1,zero_ind);


segment3 = JX_lineDensity(ind2:end)'*dx;
current3 = sum(segment3);
Itot3 = Itot1;
currentSegmentSum3 = Itot3/(N/2);
length3 = Y(end,zero_ind)-Y(ind2,zero_ind);

%Calculate optimal, even, number of compensation coil turns (assuming 1 A current):
% optimalN = round_even((current2)) %Two turns = one full turn on each side of coil.

subSegmentCurrents = [0,0,0];
i=1;
j=1;
k=1;
subStartInd = [i,j,k];
indices = {[],[],[]};
widths = {[],[],[]};
subSegmentAverages = {[],[],[]};
seedPointInd=[];
findSeedPoint=1;
for m = 1:N
    if m <= N/2
    while subSegmentCurrents(1) < abs(currentSegmentSum1)
        subSegmentCurrents(1) = abs(sum(segment1(subStartInd(1):i)));
        i=i+1;
        
        if(i>length(segment1))
            i=i-1;
        break;
        end
    end
    end
    
    while subSegmentCurrents(2) < abs(currentSegmentSum2)
        subSegmentCurrents(2) = abs(sum(segment2(subStartInd(2):j)));
        j=j+1;
        if(subSegmentCurrents(2) > abs(currentSegmentSum2/2) && findSeedPoint==1)
            seedPointInd=[seedPointInd,j];
            findSeedPoint=0;
        end
        if(j>length(segment2))
            j=j-1;
        break;
        end
    end
    findSeedPoint=1;
    
    
    if m > N/2
    while subSegmentCurrents(3) < abs(currentSegmentSum3)
        subSegmentCurrents(3) = abs(sum(segment3(subStartInd(3):k)));
        k=k+1;
        if(k>length(segment3))
            k=k-1;
        break;
        end
    end
    k=k-1;
    end
    

if m > N/2
    widths(2:3) = {[widths{2},Y(j,zero_ind)-Y(subStartInd(2),zero_ind)],[widths{3},Y(k,zero_ind)-Y(subStartInd(3),zero_ind)]};
    indices(2:3) = {[indices{2},j],[indices{3},k]};
    subSegmentAverages(2:3) = {[subSegmentAverages{2},sum(segment2(subStartInd(2):j))/length(segment2(subStartInd(2):j))], ...
    [subSegmentAverages{3},sum(segment3(subStartInd(3):k))/length(segment3(subStartInd(3):k))]};

else
    widths(1:2) = {[widths{1},Y(i,zero_ind)-Y(subStartInd(1),zero_ind)],[widths{2},Y(j,zero_ind)-Y(subStartInd(2),zero_ind)]};
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
plot(Y(:,zero_ind),(JX_lineDensity));
hold on;
scatter(Y([indices{1}],zero_ind),zeros(size(Y([indices{1}],zero_ind))));
scatter(Y([indices{2}+ind1],zero_ind),zeros(size(Y([indices{2}+ind1],zero_ind))));
scatter(Y([indices{3}+ind2],zero_ind),zeros(size(Y([indices{3}+ind2],zero_ind))));

[a1,b1]=unique(Y([indices{2}+ind1]));
b1(end)
indices{2}=indices{2}(1:b1(end));
I=2/numel(indices{2});
%% Calculate streamline seed points and the corresponding streamlines
%%%%%%%%%%%%%%%%%%%%%
% Option 1: Midway between calculated index points
index_list = [ind1,indices{2}+ind1];
seedPointIndices=round(mean([index_list(1:end-1);index_list(2:end)]));
seedPointInd(:)=0;
%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%
% Option 2: Geometrical center points (current split in half)
% seedPointIndices=[ind1+seedPointInd];
% seedPointIndices(1)=ind1+seedPointInd(1)-shiftPts1;
% seedPointIndices(2)=ind1+seedPointInd(2)-shiftPts2;
% seedPointIndices(3)=ind1+seedPointInd(3)-shiftPts3;
% seedPointIndices(4)=ind1+seedPointInd(4)-shiftPts4;
% seedPointIndices(5)=ind1+seedPointInd(5)-shiftPts5;
% seedPointIndices(6)=ind1+seedPointInd(6)-shiftPts6;
% seedPointIndices(7)=ind1+seedPointInd(7)-shiftPts7;
% seedPointIndices(8)=ind1+seedPointInd(8)-shiftPts8;
% seedPointIndices(9)=ind1+seedPointInd(9)-shiftPts9;
% seedPointIndices(10)=ind1+seedPointInd(10)-shiftPts10;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

shiftPts1=0;
shiftPts2=0;
shiftPts3=7;
shiftPts4=5;
shiftPts5=12;
shiftPts6=8;
shiftPts7=0;
shiftPts8=0;
shiftPts9=0;
shiftPts10=0;
%Manually shift slightly to get sensible streamlines:
seedPointIndices(1)=seedPointIndices(1)-shiftPts1;
seedPointIndices(2)=seedPointIndices(2)-shiftPts2;
seedPointIndices(3)=seedPointIndices(3)-shiftPts3;
seedPointIndices(4)=seedPointIndices(4)-shiftPts4;
seedPointIndices(5)=seedPointIndices(5)-shiftPts5;
seedPointIndices(6)=seedPointIndices(6)-shiftPts6;
% seedPointIndices(7)=seedPointIndices(7)-shiftPts7;
% seedPointIndices(8)=seedPointIndices(8)-shiftPts8;


[max_value,idx]=max(J_Abs(:));
[xStart,yStart,z,w]=ind2sub(size(J_Abs),idx);
streamLineList = {};
streamLineLength = 15000;
streamLineDx = 2;
cutPoints = 1;%Cutpoints * streamLineDx = length cut from streamline to facilitate interconnections

%Optional fix for streamlines that are too close to J = 0
fixedY = Y;
% fixedY(seedPointIndices(1),zero_ind) = Y(seedPointIndices(1),zero_ind) + 0.4;
% fixedY(seedPointIndices(end),zero_ind) = Y(seedPointIndices(end),zero_ind) - 0.4;

for i = 1:length(seedPointIndices)
    startx = X(seedPointIndices(i),zero_ind);
    starty = fixedY(seedPointIndices(i),zero_ind);
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
    streamLineList{i} = streamLine(1:ind(3)-cutPoints,:);
end

streamLineList

figure();
for i = 1:length(streamLineList)
    XY=streamLineList{i};
    plot(XY(:,1),XY(:,2),'linewidth',2);
    hold on;
end


streamLineLengths = [];
for i = 1:length(streamLineList)
    len = 0;
    for j = 1:length(streamLineList{i})-1
        len = len + sqrt( (streamLineList{i}(j+1,1)-streamLineList{i}(j,1))^2 + (streamLineList{i}(j+1,2)-streamLineList{i}(j,2))^2 );
    end
    streamLineLengths = [streamLineLengths,len];
end

streamLineLengths



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

maxWidth = [0.2,0.4,0.8,0.8,0.4,0.2]; %[mm], adjust widths according to each trace length
traceStandardWidth=0.5;
Npts = 101; %Should be an odd number
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
    pts = linspace(-maxWidth(i)/2,maxWidth(i)/2,Npts);
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
        if j > 1
            if j==2
            plot(streamLine(:,1),streamLine(:,2));
            end
            hold on;
            plot(integrationPathX,integrationPathY);
            axis equal;
            drawnow;
        end
        toc
        '2.'
        tic
        % 2. Interpolate Jx and Jy along the perpendicular line
        Jx_line = interp2(X,Y,JX,integrationPathX,integrationPathY);
        Jy_line = interp2(X,Y,JY,integrationPathX,integrationPathY);
        toc
        '3.'
        tic
        % 3. Take the dot product between J_line and unitVector 
        a=[Jx_line',Jy_line'];
        b=unitVector';
        dotProduct = a*b;
        dotProduct = dotProduct';
        toc
        '4. & 5.'
        tic
        % 4. and 5. Perform line integral. First starting from zero and integrating 
        % in the one direction, then starting from zero and integrating in 
        % the other direction until each integral = I/2
        integralDir1 = 0;
        integralDir2 = 0;
        startingPoint = (numel(dotProduct)+1)/2;
        m = 1;
        n = 1;
        flippedDotProduct = fliplr(dotProduct);
        while integralDir1 < currentLimit;
            if(flippedDotProduct(startingPoint+m)<0) %Stop if trying to integrate currents flowing in wrong direction
                break;
            end
            integralDir1 = sum(flippedDotProduct(startingPoint:startingPoint+m))*dl*0.001;
            m=m+1;
            if(startingPoint+m > length(dotProduct))
                m=m-1;
                break;
            elseif integralDir1 < 0
                    break;
%             elseif(i==1 || i==length(streamLine)) %optional if statement that fixes the smallest loop trace widths
%                 if integralDir1 > 0.7* I/N
%                     break;
%                 end
            end
        end
        while integralDir2 < currentLimit;
            if(flippedDotProduct(startingPoint+n)<0) %Stop if trying to integrate currents flowing in wrong direction
                break;
            end
            integralDir2 = sum(dotProduct(startingPoint:startingPoint+n))*dl*0.001;
            n=n+1;
            if(startingPoint+n > length(dotProduct))
                n=n-1;
                break;
            elseif integralDir2 < 0
                    break;
%             elseif(i==1 || i==length(streamLine)) %optional if statement that fixes the smallest loop trace widths
%                 if integralDir2 > 0.7* I/N
%                     break;
%                 end
            end
        end
        toc
        %6. Save coordinates of trace outline
        traceOutlines{i}(j,:) = [integrationPathX(startingPoint-m),integrationPathY(startingPoint-m),integrationPathX(startingPoint+n),integrationPathY(startingPoint+n)];
        
        
        
    end
    
%7. Optional: Find splitting lines for slotting conductor
% for q = 1:length(streamLine)    
%     if i<length(streamLineList)/2+1
%         x1=traceOutlines{i}(q,1);
%         x2=traceOutlines{i}(q,3);
%         y1=traceOutlines{i}(q,2);
%         y2=traceOutlines{i}(q,4);
%     else
%         x2=traceOutlines{i}(q,1);
%         x1=traceOutlines{i}(q,3);
%         y2=traceOutlines{i}(q,2);
%         y1=traceOutlines{i}(q,4);
%     end
%         localTraceWidth=sqrt((x2-x1)^2+(y2-y1)^2);
%         unitVector1 = [x2-x1,y2-y1]/localTraceWidth;
%         
%         num=1;
%         if localTraceWidth > num*splitWidthLimit
%             try
%                 splitLines{i}{num}(end+1,1:2)=[x1,y1]+unitVector1*num*traceStandardWidth;
%             catch
%                 splitLines{i}{num}=[x1,y1]+unitVector1*num*traceStandardWidth;
%             end
%         end
%         num=2;
%         if localTraceWidth > num*splitWidthLimit
%             try
%                 splitLines{i}{num}(end+1,1:2)=[x1,y1]+unitVector1*num*traceStandardWidth;
%             catch
%                 splitLines{i}{num}=[x1,y1]+unitVector1*num*traceStandardWidth;
%             end
%         end
%         num=3;
%         if localTraceWidth > num*splitWidthLimit
%             try
%                 splitLines{i}{num}(end+1,1:2)=[x1,y1]+unitVector1*num*traceStandardWidth;
%             catch
%                 splitLines{i}{num}=[x1,y1]+unitVector1*num*traceStandardWidth;
%             end
%         end
%         num=4;
%         if localTraceWidth > num*splitWidthLimit
%             try
%                 splitLines{i}{num}(end+1,1:2)=[x1,y1]+unitVector1*num*traceStandardWidth;
%             catch
%                 splitLines{i}{num}=[x1,y1]+unitVector1*num*traceStandardWidth;
%             end
%         end
%         num=5;
%         if localTraceWidth > num*splitWidthLimit
%             try
%                 splitLines{i}{num}(end+1,1:2)=[x1,y1]+unitVector1*num*traceStandardWidth;
%             catch
%                 splitLines{i}{num}=[x1,y1]+unitVector1*num*traceStandardWidth;
%             end
%         end
%         num=6;
%         if localTraceWidth > num*splitWidthLimit
%             try
%                 splitLines{i}{num}(end+1,1:2)=[x1,y1]+unitVector1*num*traceStandardWidth;
%             catch
%                 splitLines{i}{num}=[x1,y1]+unitVector1*num*traceStandardWidth;
%             end
%         end
%         num=7;
%         if localTraceWidth > num*splitWidthLimit
%             try
%                 splitLines{i}{num}(end+1,1:2)=[x1,y1]+unitVector1*num*traceStandardWidth;
%             catch
%                 splitLines{i}{num}=[x1,y1]+unitVector1*num*traceStandardWidth;
%             end
%         end
%         num=8;
%         if localTraceWidth > num*splitWidthLimit
%             try
%                 splitLines{i}{num}(end+1,1:2)=[x1,y1]+unitVector1*num*traceStandardWidth;
%             catch
%                 splitLines{i}{num}=[x1,y1]+unitVector1*num*traceStandardWidth;
%             end
%         end
%         num=9;
%         if localTraceWidth > num*splitWidthLimit
%             try
%                 splitLines{i}{num}(end+1,1:2)=[x1,y1]+unitVector1*num*traceStandardWidth;
%             catch
%                 splitLines{i}{num}=[x1,y1]+unitVector1*num*traceStandardWidth;
%             end
%         end
%         num=10;
%         if localTraceWidth > num*splitWidthLimit
%             try
%                 splitLines{i}{num}(end+1,1:2)=[x1,y1]+unitVector1*num*traceStandardWidth;
%             catch
%                 splitLines{i}{num}=[x1,y1]+unitVector1*num*traceStandardWidth;
%             end
%         end
%         num=11;
%         if localTraceWidth > num*splitWidthLimit
%             try
%                 splitLines{i}{num}(end+1,1:2)=[x1,y1]+unitVector1*num*traceStandardWidth;
%             catch
%                 splitLines{i}{num}=[x1,y1]+unitVector1*num*traceStandardWidth;
%             end
%         end
%         num=12;
%         if localTraceWidth > num*splitWidthLimit
%             try
%                 splitLines{i}{num}(end+1,1:2)=[x1,y1]+unitVector1*num*traceStandardWidth;
%             catch
%                 splitLines{i}{num}=[x1,y1]+unitVector1*num*traceStandardWidth;
%             end
%         end
%         num=13;
%         if localTraceWidth > num*splitWidthLimit
%             try
%                 splitLines{i}{num}(end+1,1:2)=[x1,y1]+unitVector1*num*traceStandardWidth;
%             catch
%                 splitLines{i}{num}=[x1,y1]+unitVector1*num*traceStandardWidth;
%             end
%         end
%         num=14;
%         if localTraceWidth > num*splitWidthLimit
%             try
%                 splitLines{i}{num}(end+1,1:2)=[x1,y1]+unitVector1*num*traceStandardWidth;
%             catch
%                 splitLines{i}{num}=[x1,y1]+unitVector1*num*traceStandardWidth;
%             end
%         end
%         num=15;
%         if localTraceWidth > num*splitWidthLimit
%             try
%                 splitLines{i}{num}(end+1,1:2)=[x1,y1]+unitVector1*num*traceStandardWidth;
%             catch
%                 splitLines{i}{num}=[x1,y1]+unitVector1*num*traceStandardWidth;
%             end
%         end
%         num=16;
%         if localTraceWidth > num*splitWidthLimit
%             try
%                 splitLines{i}{num}(end+1,1:2)=[x1,y1]+unitVector1*num*traceStandardWidth;
%             catch
%                 splitLines{i}{num}=[x1,y1]+unitVector1*num*traceStandardWidth;
%             end
%         end
%         num=17;
%         if localTraceWidth > num*splitWidthLimit
%             try
%                 splitLines{i}{num}(end+1,1:2)=[x1,y1]+unitVector1*num*traceStandardWidth;
%             catch
%                 splitLines{i}{num}=[x1,y1]+unitVector1*num*traceStandardWidth;
%             end
%         end
%         num=18;
%         if localTraceWidth > num*splitWidthLimit
%             try
%                 splitLines{i}{num}(end+1,1:2)=[x1,y1]+unitVector1*num*traceStandardWidth;
%             catch
%                 splitLines{i}{num}=[x1,y1]+unitVector1*num*traceStandardWidth;
%             end
%         end
%         
% end
% 
% try
%     for ii=1:num
%     plot(splitLines{i}{ii}(:,1),splitLines{i}{ii}(:,2),'r','linewidth',1);
%     end
% catch
% end
% End 7.


    color=[rand,rand,rand];
    plot(traceOutlines{i}(:,1),traceOutlines{i}(:,2),'color','b','linewidth',1);
    hold on;
    plot(traceOutlines{i}(:,3),traceOutlines{i}(:,4),'color','b','linewidth',1);
    XY=streamLineList{i};
    plot(XY(:,1),XY(:,2),'black','linewidth',1);
    drawnow;
end


%% OPTIONAL: Split thick conductors into thinner parallel conductors

splitWidthLimit=0.8;
offset1=0.0; %push divisions forward to facilitate connections
split11=[];
split12=[];
split2=[];
split2_11=[];
split2_12=[];
split3=[];
split4=[];
split5=[];
split6=[];

conductorNumber=1;

for i = conductorNumber%length(traceOutlines) %Iterate only through a single conductor at a time
    x1=traceOutlines{i}(:,1);
    y1=traceOutlines{i}(:,2);
    x2=traceOutlines{i}(:,3);
    y2=traceOutlines{i}(:,4);
    for j=1:length(x1)
        
        Dx=(x2(j)-x1(j));
        Dy=(y2(j)-y1(j));
        localTraceWidth=sqrt(Dx^2+Dy^2)
        if localTraceWidth > 14*splitWidthLimit
            
        elseif localTraceWidth > 13*splitWidthLimit
            
        elseif localTraceWidth > (11+offset1)*splitWidthLimit && localTraceWidth < 13*splitWidthLimit
            lim=12;
            split6(end+1,1:22)=[...
                x1(j)+(x2(j)-x1(j))/lim,y1(j)+(y2(j)-y1(j))/lim,...
                x1(j)+2*Dx/lim,y1(j)+2*Dy/lim,...
                x1(j)+3*Dx/lim,y1(j)+3*Dy/lim,...
                x1(j)+4*Dx/lim,y1(j)+4*Dy/lim,...
                x1(j)+5*Dx/lim,y1(j)+5*Dy/lim,...
                x1(j)+6*Dx/lim,y1(j)+6*Dy/lim,...
                x1(j)+7*Dx/lim,y1(j)+7*Dy/lim,...
                x1(j)+8*Dx/lim,y1(j)+8*Dy/lim,...
                x1(j)+9*Dx/lim,y1(j)+9*Dy/lim,...
                x1(j)+10*Dx/lim,y1(j)+10*Dy/lim,...
                x1(j)+11*Dx/lim,y1(j)+11*Dy/lim];  
        elseif localTraceWidth > (9+offset1)*splitWidthLimit && localTraceWidth < 11*splitWidthLimit
            lim=10;
            split5(end+1,1:18)=[...
                x1(j)+(x2(j)-x1(j))/lim,y1(j)+(y2(j)-y1(j))/lim,...
                x1(j)+2*Dx/lim,y1(j)+2*Dy/lim,...
                x1(j)+3*Dx/lim,y1(j)+3*Dy/lim,...
                x1(j)+4*Dx/lim,y1(j)+4*Dy/lim,...
                x1(j)+5*Dx/lim,y1(j)+5*Dy/lim,...
                x1(j)+6*Dx/lim,y1(j)+6*Dy/lim,...
                x1(j)+7*Dx/lim,y1(j)+7*Dy/lim,...
                x1(j)+8*Dx/lim,y1(j)+8*Dy/lim,...
                x1(j)+9*Dx/lim,y1(j)+9*Dy/lim];  
        elseif localTraceWidth > (7+offset1)*splitWidthLimit && localTraceWidth < 9*splitWidthLimit
            lim=8;
            split4(end+1,1:14)=[...
                x1(j)+Dx/lim,y1(j)+Dy/lim,...
                x1(j)+2*Dx/lim,y1(j)+2*Dy/lim,...
                x1(j)+3*Dx/lim,y1(j)+3*Dy/lim,...
                x1(j)+4*Dx/lim,y1(j)+4*Dy/lim,...
                x1(j)+5*Dx/lim,y1(j)+5*Dy/lim,...
                x1(j)+6*Dx/lim,y1(j)+6*Dy/lim,...
                x1(j)+7*Dx/lim,y1(j)+7*Dy/lim];   
        elseif localTraceWidth > (5+offset1)*splitWidthLimit && localTraceWidth < 7*splitWidthLimit
            split3(end+1,1:10)=[x1(j)+Dx/6,y1(j)+Dy/6,...
                x1(j)+2*Dx/6,y1(j)+2*Dy/6,...
                x1(j)+3*Dx/6,y1(j)+3*Dy/6,...
                x1(j)+4*Dx/6,y1(j)+4*Dy/6,...
                x1(j)+5*Dx/6,y1(j)+5*Dy/6];
        elseif localTraceWidth > (3+offset1)*splitWidthLimit && localTraceWidth < 5*splitWidthLimit
            split2(end+1,1:6)=[x1(j)+Dx/4,y1(j)+Dy/4,...
                x1(j)+2*Dx/4,y1(j)+2*Dy/4,...
                x1(j)+3*Dx/4,y1(j)+3*Dy/4];
            split2_11(end+1,1:2)=[split2(end,5)+Dx*(localTraceWidth-splitWidthLimit)/12/localTraceWidth,split2(end,6)+Dy*(localTraceWidth-splitWidthLimit)/12/localTraceWidth];
            split2_12(end+1,1:2)=[split2(end,5)-Dx*(localTraceWidth-splitWidthLimit)/12/localTraceWidth,split2(end,6)-Dy*(localTraceWidth-splitWidthLimit)/12/localTraceWidth];
        elseif localTraceWidth > 1*splitWidthLimit && localTraceWidth < 3*splitWidthLimit
            split11(end+1,1:2)=[x1(j)+Dx*splitWidthLimit/(2*localTraceWidth),y1(j)+Dy*splitWidthLimit/(2*localTraceWidth)];
            split12(end+1,1:2)=[x2(j)-Dx*splitWidthLimit/(2*localTraceWidth),y2(j)-Dy*splitWidthLimit/(2*localTraceWidth)];
        end
    end
end

figure();
plot(traceOutlines{conductorNumber}(:,1),traceOutlines{conductorNumber}(:,2));
hold on;
plot(traceOutlines{conductorNumber}(:,3),traceOutlines{conductorNumber}(:,4));
plot(split11(:,1),split11(:,2));
plot(split12(:,1),split12(:,2));

plot(split2(:,1),split2(:,2));
plot(split2(:,3),split2(:,4));
plot(split2(:,5),split2(:,6));

plot(split2_11(:,1),split2_11(:,2),'r');
plot(split2_12(:,1),split2_12(:,2),'r');

plot(split3(:,1),split3(:,2));
plot(split3(:,3),split3(:,4));
plot(split3(:,5),split3(:,6));
plot(split3(:,7),split3(:,8));
plot(split3(:,9),split3(:,10));

plot(split4(:,1),split4(:,2));
plot(split4(:,3),split4(:,4));
plot(split4(:,5),split4(:,6));
plot(split4(:,7),split4(:,8));
plot(split4(:,9),split4(:,10));
plot(split4(:,11),split4(:,12));
plot(split4(:,13),split4(:,14));

plot(split5(:,1),split5(:,2));
plot(split5(:,3),split5(:,4));
plot(split5(:,5),split5(:,6));
plot(split5(:,7),split5(:,8));
plot(split5(:,9),split5(:,10));
plot(split5(:,11),split5(:,12));
plot(split5(:,13),split5(:,14));
plot(split5(:,15),split5(:,16));
plot(split5(:,17),split5(:,18));

plot(split6(:,1),split6(:,2));
plot(split6(:,3),split6(:,4));
plot(split6(:,5),split6(:,6));
plot(split6(:,7),split6(:,8));
plot(split6(:,9),split6(:,10));
plot(split6(:,11),split6(:,12));
plot(split6(:,13),split6(:,14));
plot(split6(:,15),split6(:,16));
plot(split6(:,17),split6(:,18));
plot(split6(:,19),split6(:,20));
plot(split6(:,21),split6(:,22));

%% Output individual polygons

output1=[];
output2=[];
output3=[];
output4=[];
output5=[];
output6=[];
finalpoint=[];
output1 = [output1;traceOutlines{1}(:,1:2)];
output1 = [output1;fliplr(traceOutlines{1}(:,3:4)')'];
output1 = [output1;output1(1,:)];
output2 = [output2;traceOutlines{2}(:,1:2)];
output2 = [output2;fliplr(traceOutlines{2}(:,3:4)')'];
output2 = [output2;output2(1,:)];
output3 = [output3;traceOutlines{3}(:,1:2)];
output3 = [output3;fliplr(traceOutlines{3}(:,3:4)')'];
output3 = [output3;output3(1,:)];
output4 = [output4;traceOutlines{4}(:,1:2)];
output4 = [output4;fliplr(traceOutlines{4}(:,3:4)')'];
output4 = [output4;output4(1,:)];
output5 = [output5;traceOutlines{5}(:,1:2)];
output5 = [output5;fliplr(traceOutlines{5}(:,3:4)')'];
output5 = [output5;output5(1,:)];
output6 = [output6;traceOutlines{6}(:,1:2)];
output6 = [output6;fliplr(traceOutlines{6}(:,3:4)')'];
output6 = [output6;output6(1,:)];

%% Concatenate all data to be exported as a single polygon

output=[];
finalpoint=[];
for i = 1:length(traceOutlines)
   
    output = [output;traceOutlines{i}(:,1:2)];
    
end
for i = length(traceOutlines):-1:1

    output = [output;fliplr(traceOutlines{i}(:,3:4)')'];
    
end
output = [output;finalpoint];
output = [output;output(1,:)];
figure();
plot(output(:,1),output(:,2))

%% Change directions of connections if necessary
output=[];
finalpoint=[];
for i = 1:length(traceOutlines)
   
    output = [output;fliplr(traceOutlines{i}(:,1:2)')'];
    
end
for i = length(traceOutlines):-1:1

    output = [output;traceOutlines{i}(:,3:4)];
    
end
output = [output;finalpoint];
output = [output;output(1,:)];
figure();
plot(output(:,1),output(:,2))

%% Only two loops
output=[];
finalpoint=[];
for i = 2:2
   
    output = [output;fliplr(traceOutlines{i}(:,1:2)')'];
    
end

for i = 2:-1:2

    output = [output;traceOutlines{i}(:,3:4)];
    
end
