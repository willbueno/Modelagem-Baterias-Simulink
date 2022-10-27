% Inicializa a estrutura dos dados para o simulink
PulseData = struct('voltage',[],'current',[],'time',[],'temperature',[]);

% Inicializa as variáveis para leitura dos dados experimentais de tensão
filename = 'data_voltage.csv';
delimiter = ';';
formatSpec = '%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter',...
    delimiter, 'TextType', 'string', 'EmptyValue', NaN,...
    'ReturnOnError', false);
fclose(fileID);

% Armazena os valores na estrutura
PulseData.voltage = [dataArray{1:end-1}];

% Deletar variáveis temporárias
clearvars filename delimiter formatSpec fileID dataArray ans;

% Inicializa as variáveis para leitura dos dados experimentais de corrente
filename = 'data_current.csv';
delimiter = ';';
formatSpec = '%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter',...
    delimiter, 'TextType', 'string', 'EmptyValue', NaN,...
    'ReturnOnError', false);
fclose(fileID);

% Armazena os valores na estrutura
PulseData.current = [dataArray{1:end-1}];

% Deletar variáveis temporárias
clearvars filename delimiter formatSpec fileID dataArray ans;

% Inicializa as variáveis para leitura dos dados experimentais de tempo
filename = 'data_timer.csv';
delimiter = ';';
formatSpec = '%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter',...
    delimiter, 'TextType', 'string', 'EmptyValue', NaN,...
    'ReturnOnError', false);
fclose(fileID);

% Armazena os valores na estrutura
PulseData.time = [dataArray{1:end-1}];

% Deletar variáveis temporárias
clearvars filename delimiter formatSpec fileID dataArray ans;

% Temperatura
tam = size(PulseData.time);
PulseData.temperature = 25*ones(tam-1);
clear tam

% Deletando o último valor dos vetores
PulseData.time = PulseData.time(1:end-1);
PulseData.voltage = PulseData.voltage(1:end-1);
PulseData.current = (-1)*PulseData.current(1:end-1);

% Numero de pulsos
pulsos = 100;

% Corrente do modelo
current = -1.5;

% SOC Lookup Table breakpoints
SOC_LUT = (0:0.01:1)';

% Capacidade da bateria
% Medição através da contagem de coulomb da curva de descarga
Capacity = 2; % Ah

% Capacidade drenada inicial
Qe_init = 0; %Ampere*hours

% Inicialização dos parâmetros
% Em OCV vs SOC
Em = ones(size(SOC_LUT)); % Volts

% R0 vs SOC
R0 = ones(size(SOC_LUT)); % Ohms

% R1 vs SOC
R1 = ones(size(SOC_LUT)); % Ohms

% C1 vs SOC
C1 = ones(size(SOC_LUT)); % Farad

% Encontrar os pontos de interesse
ntotal = size(PulseData.time);
ptime = zeros(2*pulsos+1,1);
pvolt = zeros(2*pulsos+1,1);
pPosition = zeros(2*pulsos+1,1);
ovolt = zeros(2*pulsos+1,1);
otime = zeros(2*pulsos+1,1);
oPosition = zeros(2*pulsos+1,1);
aux = 1;

for k=2:ntotal
    
    if PulseData.current(k-1) >= (current+0.1) && ...
            PulseData.current(k) < (current+0.1)
        
        for l=k-10:k
            
            if (PulseData.voltage(l)-PulseData.voltage(l+1)) > 0.1
                
                pvolt(aux) = PulseData.voltage(l);
                ptime(aux) = PulseData.time(l);
                pPosition(aux) = l;
                
                ovolt(aux) = PulseData.voltage(l+3);
                otime(aux) = PulseData.time(l+3);
                oPosition(aux) = l+3;
                
                aux = aux+1;
                
            end
            
        end
        
    end
    
    if PulseData.current(k-1) <= (current+0.1) && ...
            PulseData.current(k) > (current+0.1)
        
        for m=k-10:k
            
            if (PulseData.voltage(m+1)-PulseData.voltage(m)) > 0.1
                
                pvolt(aux) = PulseData.voltage(m);
                ptime(aux) = PulseData.time(m);
                pPosition(aux) = m;
                
                ovolt(aux) = PulseData.voltage(m+3);
                otime(aux) = PulseData.time(m+3);
                oPosition(aux) = m+3;
                
                aux = aux+1;
                
            end
            
        end
        
    end
    
end

% Vetores de pontos de interesse
pvolt(end) = PulseData.voltage(end);
ptime(end) = PulseData.time(end);
pPosition(end) = PulseData.time(end);
ovolt(end) = ovolt(end-1);
otime(end) = pvolt(end-1);
oPosition(end) = pPosition(end-1);

% Cálculo dos parâmetros
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
                
                % T = R*C -> C = T/R
                C1aux = (PulseData.time(l) - otime(k-1))/R1(auxR1);
                break
                
            end
            
        end
        
        for l=oPosition(k):ntotal
            
            if PulseData.voltage(l) >= C1aux2
                
                % T = R*C -> C = T/R
                C1aux2 = (PulseData.time(l) - otime(k))/R1(auxR1);
                break
                
            end
            
        end
        
        C1(auxR1) = (C1aux+C1aux2)/2;
        
        auxR1 = auxR1-1;
        
    end
    
end

% Estado de Carga 0%
R1(1) = R1(2);
R0(1) = R0(2);
C1(1) = C1(2);

% Limpar variáveis auxialiares
clear C1aux C1aux2 aux auxR1 current k l m
clear ntotal ntotal2 pulsos otime ovolt
clear ptime pvolt oPosition pPosition