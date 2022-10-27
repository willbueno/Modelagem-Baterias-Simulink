close all
clear all
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
Capacity = 4.5; %Ampere*hours

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
pulsos = 30;
ptime = zeros(2*pulsos+1,1);
pvolt = zeros(2*pulsos+1,1);
pPosition = zeros(2*pulsos+1,1);
ovolt = zeros(2*pulsos+1,1);
otime = zeros(2*pulsos+1,1);
oPosition = zeros(2*pulsos+1,1);
current = -2;
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

for k=1:72
    Em(k) = Em(72);
    C1(k) = C1(72);
    R1(k) = R1(72);
    R0(k) = R0(72);
end

%%
figure(1)
subplot(2,1,1)
plot(PulseData.time, PulseData.voltage,'b')
legend({'Tensão'},'FontSize',18)
set(gca,'FontSize',16)
ylabel('Tensão (V)','FontSize',18)
xlabel('Tempo (s)','FontSize',18)
axis([0 7895 11.8 12.8])
grid
subplot(2,1,2)
plot(PulseData.time, PulseData.current,'b')
grid
legend({'Corrente'},'FontSize',18)
set(gca,'FontSize',16)
ylabel('Corrente (A)','FontSize',18)
xlabel('Tempo (s)','FontSize',18)
axis([0 7895 -2.5 0.5])

figure(2)
subplot(2,1,1)
plot(PulseData.time, PulseData.voltage,'b')
hold on
plot(ptime,pvolt,'*r')
hold on
plot(otime,ovolt,'ok')
grid
axis([0 7895 11.8 12.8])
subplot(2,1,2)
plot(PulseData.time, PulseData.current,'b')
grid
axis([0 7895 -2.5 0.5])

figure(3)
plot(100*SOC_LUT, Em, 'b')
grid
legend({'Em'},'FontSize',22)
set(gca,'FontSize',20)
ylabel('Tensão (V)','FontSize',22)
xlabel('Estado de carga (%)','FontSize',22)
axis([70 100 12.3 12.8])

figure(4)
grid
plot(100*SOC_LUT, R0, 'b')
grid
legend({'R0'},'FontSize',22)
set(gca,'FontSize',20)
ylabel('Resistência (\Omega)','FontSize',22)
xlabel('Estado de carga (%)','FontSize',22)
axis([70 100 0.13 0.15])

figure(5)
plot(100*SOC_LUT, R1, 'b')
grid
legend({'R1'},'FontSize',22)
set(gca,'FontSize',20)
ylabel('Resistência (\Omega)','FontSize',22)
xlabel('Estado de carga (%)','FontSize',22)
axis([70 100 0.1 0.2])

figure(6)
plot(100*SOC_LUT, C1, 'b')
grid
legend({'C1'},'FontSize',22)
set(gca,'FontSize',20)
ylabel('Capacitância (F)','FontSize',22)
xlabel('Estado de carga (%)','FontSize',22)
axis([70 100 20 70])

figure(11)
plot(100*SOC_LUT, C1.*R1, 'b')
grid
legend({'\tau'},'FontSize',22)
set(gca,'FontSize',20)
ylabel('Constante de Tempo \tau (s)','FontSize',22)
xlabel('Estado de carga (%)','FontSize',22)
axis([70 100 0 10])

%% Limpar variáveis auxialiares
clear C1aux C1aux2 aux auxR1 current k l m
clear ntotal ntotal2 pulsos otime ovolt
clear ptime pvolt oPosition pPosition