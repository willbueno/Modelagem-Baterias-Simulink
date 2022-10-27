close all
clear
clc
PulseData = struct('voltage',[],'current',[],'time',[],'temperature',[]);

%% Initialize variables.
filename = 'data_voltage.csv';
delimiter = ';';
formatSpec = '%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter',...
    delimiter, 'TextType', 'string', 'EmptyValue', NaN,  'ReturnOnError', false);
fclose(fileID);

%% Create output variable
PulseData.voltage = [dataArray{1:end-1}];

%% Clear temporary variables
clearvars filename delimiter formatSpec fileID dataArray ans;


%% Initialize variables.
filename = 'data_current.csv';
delimiter = ';';
formatSpec = '%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter',...
    delimiter, 'TextType', 'string', 'EmptyValue', NaN,  'ReturnOnError', false);
fclose(fileID);

%% Create output variable
PulseData.current = [dataArray{1:end-1}];

%% Clear temporary variables
clearvars filename delimiter formatSpec fileID dataArray ans;

%% Initialize variables.
filename = 'data_timer.csv';
delimiter = ';';
formatSpec = '%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter',...
    delimiter, 'TextType', 'string', 'EmptyValue', NaN,  'ReturnOnError', false);
fclose(fileID);

%% Create output variable
PulseData.time = [dataArray{1:end-1}];

%% Clear temporary variables
clearvars filename delimiter formatSpec fileID dataArray ans;

%% Temperatura
tam = size(PulseData.time);
PulseData.temperature = 25*ones(tam);
clear tam

%% Ajustes
PulseData.time = PulseData.time(1:end-1);
PulseData.voltage = PulseData.voltage(1:end-1);
PulseData.current = (-1)*PulseData.current(1:end-1);

%% Chosen Values
% SOC Lookup Table breakpoints
SOC_LUT = (0:0.01:1)';

%% Known Values
% Battery capacity
% Measured by coulomb counting the discharge curve
Capacity = 2; %Ampere*hours

% Charge deficit at start of data set
% Assumption based on preparation of test
Qe_init = 0; %Ampere*hours

%% Estimated Parameters - Initial starting points before estimation
% Em open-circuit voltage vs SOC
Em = ones(size(SOC_LUT)); %Volts

% R0 resistance vs SOC
R0 = ones(size(SOC_LUT)); %Ohms

% R1 Resistance vs SOC
R1 = ones(size(SOC_LUT)); %Ohms

% C1 Capacitance vs SOC
C1 = ones(size(SOC_LUT)); %Farad

%% Encontrar os pontos de interesse
ntotal = size(PulseData.time);
pulsos = 100;
ptime = zeros(2*pulsos+1,1);
pvolt = zeros(2*pulsos+1,1);
pPosition = zeros(2*pulsos+1,1);
ovolt = zeros(2*pulsos+1,1);
otime = zeros(2*pulsos+1,1);
oPosition = zeros(2*pulsos+1,1);
current = -1.5;
aux = 1;

for k=2:ntotal
    if PulseData.current(k-1) >=  (current+0.1) && PulseData.current(k) < (current+0.1)
        for l=k-10:k
            if (PulseData.voltage(l)-PulseData.voltage(l+1)) > 0.1
                pvolt(aux) = PulseData.voltage(l);
                ptime(aux) = PulseData.time(l);
                pPosition(aux) = l;
                
                ovolt(aux) = PulseData.voltage(l+3); % Arrumar para pegar um valor na transição por codigo
                otime(aux) = PulseData.time(l+3); % Arrumar para pegar um valor na transição por codigo
                oPosition(aux) = l+3;
                
                aux = aux+1;
            end
        end
    end
    
    if PulseData.current(k-1) <=  (current+0.1) && PulseData.current(k) > (current+0.1)
        for m=k-10:k
            if (PulseData.voltage(m+1)-PulseData.voltage(m)) > 0.1
                pvolt(aux) = PulseData.voltage(m);
                ptime(aux) = PulseData.time(m);
                pPosition(aux) = m;
                
                ovolt(aux) = PulseData.voltage(m+3); % Arrumar para pegar um valor na transição por codigo
                otime(aux) = PulseData.time(m+3); % Arrumar para pegar um valor na transição por codigo
                oPosition(aux) = m+3;
                
                aux = aux+1;
            end
        end
    end
end

pvolt(end) = PulseData.voltage(end);
ptime(end) = PulseData.time(end);
pPosition(end) = PulseData.time(end);
ovolt(end) = ovolt(end-1);
otime(end) = pvolt(end-1);
oPosition(end) = pPosition(end-1);

%% Calculo dos parametros
ntotal2 = size(pvolt);
aux = 101;
auxR1 = 101;

for k=1:ntotal2
    if rem(k,2) == 1
        Em(aux) = pvolt(k);
        R0(aux) = (pvolt(k)-ovolt(k))/((-1)*current);
        aux = aux-1;
    end
    if rem(k,2) == 0
        R1(auxR1) = (ovolt(k-1)-pvolt(k))/((-1)*current);
        C1aux = ovolt(k-1)-(0.632*(ovolt(k-1)-pvolt(k)));
        C1aux2 = ovolt(k)+(0.632*(pvolt(k+1)-ovolt(k)));
        
        for l=oPosition(k-1):ntotal
            if PulseData.voltage(l) <= C1aux
                C1aux = (PulseData.time(l) - otime(k-1))/R1(auxR1); % T = R*C -> C = T/R
                break
            end
        end
        
        for l=oPosition(k):ntotal
            if PulseData.voltage(l) >= C1aux2
                C1aux2 = (PulseData.time(l) - otime(k))/R1(auxR1); % T = R*C -> C = T/R
                break
            end
        end

        C1(auxR1) = (C1aux+C1aux2)/2;
        
        auxR1 = auxR1-1;
    end
end

R1(1) = R1(2);
R0(1) = R0(2);
C1(1) = C1(2);

%% Plot variaveis do modelo
figure(1)
subplot(2,1,1)
plot(PulseData.time, PulseData.voltage,'b')
legend('Tensão')
ylabel('Tensão (V)')
xlabel('Tempo (s)')
axis([0 22300 10 12.5])
grid
subplot(2,1,2)
plot(PulseData.time, PulseData.current,'b')
grid
legend('Corrente')
ylabel('Corrente (A)')
xlabel('Tempo (s)')
axis([0 22300 -2 0.5])

figure(2)
subplot(2,1,1)
plot(PulseData.time, PulseData.voltage,'b')
hold on
plot(ptime,pvolt,'*r')
hold on
plot(otime,ovolt,'ok')
grid
subplot(2,1,2)
plot(PulseData.time, PulseData.current,'b')
grid

figure(3)
plot(100*SOC_LUT, Em, 'b')
grid
legend('Em')
ylabel('Tensão (V)')
xlabel('Estado de carga (%)')

figure(4)
grid
plot(100*SOC_LUT, R0, 'b')
grid
legend('R0')
ylabel('Resistência (?)')
xlabel('Estado de carga (%)')

figure(5)
plot(100*SOC_LUT, R1, 'b')
grid
legend('R1')
ylabel('Resistência (?)')
xlabel('Estado de carga (%)')

figure(6)
plot(100*SOC_LUT, C1, 'b')
grid
legend('C1')
ylabel('Capacitância (F)')
xlabel('Estado de carga (%)')

figure(7)
plot(PulseData.time, PulseData.voltage,'b')
hold on
plot(ptime,pvolt,'*r')
hold on
plot(otime,ovolt,'ok')
grid
axis([0 233.1 11.5 13])
annotation('textarrow',[0.17 0.15],[0.7 0.65],'String','E_m @ SOC_a')
annotation('textarrow',[0.88 0.905],[0.66 0.625],'String','E_m @ SOC_b')
annotation('line',[0.157 0.2], [0.375 0.375])
annotation('line',[0.157 0.2], [0.645 0.645])
annotation('doublearrow',[0.1785 0.1785],[0.375 0.645])
text(20, 12.23, {'Queda de tensão', 'R_0 @ SOC_a'})
annotation('textarrow',[0.22 0.16],[0.2 0.33],'String','R-C @ SOC_a')
annotation('arrow',[0.22 0.21],[0.2 0.31])
annotation('arrow',[0.22 0.28],[0.2 0.29])
annotation('line',[0.317 0.36], [0.295 0.295])
annotation('line',[0.317 0.36], [0.565 0.565])
annotation('doublearrow',[0.3385 0.3385],[0.295 0.565])
text(70, 12.1, {'Queda de tensão', 'R_0 @ SOC_b'})
annotation('textarrow',[0.4 0.4],[0.68 0.61],'String','R-C @ SOC_b')
annotation('arrow',[0.4 0.32],[0.68 0.595])
annotation('arrow',[0.4 0.55],[0.68 0.625])
legend('Tensão')
ylabel('Tensão (V)')
xlabel('Tempo (s)')

figure(8)
plot(PulseData.time, PulseData.voltage,'b')
hold on
plot(ptime,pvolt,'*r')
hold on
plot(otime,ovolt,'ok')
grid
axis([0 233.1 11.5 13])
annotation('textarrow',[0.17 0.15],[0.7 0.65],'String','A')
annotation('textarrow',[0.185 0.157], [0.415 0.375], 'String', 'B')
annotation('textarrow',[0.34 0.312], [0.335 0.298], 'String', 'C')
annotation('textarrow',[0.28 0.3],[0.6 0.565],'String','D')
annotation('textarrow',[0.88 0.90],[0.67 0.635],'String','E')
legend('Tensão')
ylabel('Tensão (V)')
xlabel('Tempo (s)')

figure(9)
subplot(2,1,1)
plot(PulseData.time, PulseData.voltage,'b')
hold on
plot(ptime,pvolt,'*r')
hold on
plot(otime,ovolt,'ok')
grid
axis([0 461.05 11.5 13])
legend('Tensão', 'Derivada Infinita', 'Início da Dinâmica')
ylabel('Tensão (V)')
xlabel('Tempo (s)')

subplot(2,1,2)
plot(PulseData.time, PulseData.current,'b')
grid
axis([0 461.05 -2 0.5])
legend('Corrente')
ylabel('Corrente (V)')
xlabel('Tempo (s)')

%% Limpar variaveis
clear C1aux C1aux2 aux auxR1 current k l m ntotal ntotal2 pulsos otime ovolt ptime pvolt oPosition pPosition