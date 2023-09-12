

clear all
close all 
clc

filename="C:\Users\objoerkqvist\OneDrive - ETH Zurich\Documents\MATLAB\Coil design\harvesting_coil_v7_impedance_96PARAM_SWEEP.txt";
fid = fopen(filename);
tline = fgetl(fid);
runCount = 0;
params = zeros(96,3);
effective_cross_section = zeros(96,1);
emf = zeros(96,1001);
R = zeros(96,1001);
ii=1;
freq=zeros(1001,1);
B1 = 1e-6;
while ischar(tline)
    if(tline(1:11) == '#Parameters')
        runCount = runCount + 1
        params(runCount,:)=GetParams(tline);
        effective_cross_section(runCount) = params(runCount,2)*1e-3*params(runCount,3)*0.029;
        tline = fgetl(fid);
        tline = fgetl(fid);
        ii=1;
    else
        if(ii < 1002)
            aa=split(tline);
            if(runCount == 1)
                freq(ii)=str2num(aa{1});
            end
            R(runCount,ii)=str2num(aa{2});
            emf(runCount,ii) = effective_cross_section(runCount) .* freq(ii) * 1e9 * 2 * pi * B1;
            
        else

        end
        ii=ii+1; 

    end    
    tline = fgetl(fid);
end



fclose(fid);
%%

varNames = ["Wire gauge = ","Coil height = ","Number of turns = "];
user_input = [0.5,20,NaN];
sweepVarInd = find(isnan(user_input) == 1);
plotDataIndices = FindIndices(user_input,params); %Choose data to plot ([wg,H,N]): set to NaN to sweep that variable
figure;
for i = 1:numel(plotDataIndices)
    plot(freq,emf(plotDataIndices(i),:).^2 ./ R(plotDataIndices(i),:));
    hold on;
end
title(['Wire gauge = ', num2str(user_input(1)), ', Coil height = ', num2str(user_input(2)), ', Nturns = ', num2str(user_input(3))]);
legend([strcat(repmat(varNames(sweepVarInd),[numel(plotDataIndices),1]),  num2str(params(plotDataIndices(:),sweepVarInd)))])
ylabel('Power [a.u.]');
xlabel('Frequency [GHz]');
grid on;

function ind = FindIndices(input,params)
    ind = [];
    for i = 1:numel(params(:,1))
        if(params(i,1) == input(1) || isnan(input(1)))
            if(params(i,2) == input(2) || isnan(input(2)))
                if(params(i,3) == input(3) || isnan(input(3)))
                    ind = [ind,i];
                end
            end
        end
    end
end

function params = GetParams(line)

    wg_ind1 = strfind(line,'mainCoil_WG') + 12;
    wg_ind2 = strfind(line,'; mainCoil_winding_gap')-1;
    wg = str2num(line(wg_ind1:wg_ind2));

    H_ind1 = strfind(line,'mainCoil_H=') + 11;
    H_ind2 = strfind(line,'; mainCoil_L')-1;
    H = str2num(line(H_ind1:H_ind2));

    N_ind1 = strfind(line,'Nturns=') + 7;
    N_ind2 = strfind(line,'; blend_radius')-1;
    N = str2num(line(N_ind1:N_ind2));

    params=[wg,H,N];

end