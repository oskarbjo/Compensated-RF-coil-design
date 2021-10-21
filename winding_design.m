
dataPath = "C:\Users\objoerkqvist\Documents\CST\Harvesting coil designs\vertical shielded coil\vertical shielded coil multiturn\Export\isolines_0p25mmWG.txt"
rawData=fopen(dataPath)
tline = fgetl(rawData);

isoLines = {}
isoline=[]
ind=-2
lineNr = 1
while ischar(tline)
    lineNr = lineNr+1;
    rowData=str2double(strsplit(tline));
    try
        xy = rowData(2:3);
        if isnan(xy(1)) || isnan(xy(2))
            error('Not a number');
        end
        isoline=[isoline;double(rowData(2:3))];
    catch
        sprintf('new isoline found %.15g',ind)
        if ind > 0
            isoLines{ind} = isoline;
            isoline = [];
        end
        ind=ind+1;
    end
   tline=fgetl(rawData);
end
%%

figure();

for j = 1:length(isoLines)
    a=cell2mat(isoLines(j));
    plot(a(:,1),a(:,2),'linewidth',4);
    hold on;
end
    
    
%%

%[2,5,9,10,11,12,13,14,15,16]
figure()
% for j=10:length(isoLines)
for j = [15,11,9,22,10,2]
    a=cell2mat(isoLines(j));
    plot(a(:,1),a(:,2),'linewidth',4);

    hold on;

    clear a
end
legend('1','2','3','4','5','6','7','8','9','10')

%%
indices = [15,11,9,22,10,2];
a1 = cell2mat(isoLines(indices(1)));
a2 = cell2mat(isoLines(indices(2)));
a3 = cell2mat(isoLines(indices(3)));
a4 = cell2mat(isoLines(indices(4)));
a5 = cell2mat(isoLines(indices(5)));
a6 = cell2mat(isoLines(indices(6)));
a1 = a1(90:205,:);
a2 = a2(85:282,:);
a3 = a3(42:210,:);
a4 = a4(135:338,:);
a5 = a5(71:329,:);
a6 = a6(130:376,:);

%shift y positions slightly:
a2(:,2) = a2(:,2)-0.1;
a4(:,2) = a4(:,2)-0.07;
a6(:,2) = a6(:,2)+0.03;

z1 = -ones(length(a1),1);
z2 = -ones(length(a2),1);
z3 = -ones(length(a3),1);
z4 = -ones(length(a4),1);
z5 = -ones(length(a5),1);
z6 = -ones(length(a6),1);

figure()
plot(a1(:,1),a1(:,2),'linewidth',4);
hold on
plot(a2(:,1),a2(:,2),'linewidth',4);
plot(a3(:,1),a3(:,2),'linewidth',4);
plot(a4(:,1),a4(:,2),'linewidth',4);
plot(a5(:,1),a5(:,2),'linewidth',4);
plot(a6(:,1),a6(:,2),'linewidth',4);


full = [a6;a5;a4;a3;a2;a1];

figure()
plot(full(:,1),full(:,2));


a1 = [a1,z1];
a2 = [a2,z2];
a3 = [a3,z3];
a4 = [a4,z4];
a5 = [a5,z5];
a6 = [a6,z6];

full = [a6;a5;a4;a3;a2;a1];