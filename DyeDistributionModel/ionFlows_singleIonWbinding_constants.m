%% Constants and parameters of ion flow model 
           
%% Constants of use 

%Constants of use
R = 8.3144598*1e-6;     %Universal Gas Constant expressed here in units [J/umol.K = L.Pa/mol.K]. Note that R=Na*k, where k is Boltzmann coeff.
T = 298.15;             %temperature in units [K] PS: 298.15 is equivalent to 25C
F = 96485.332*1e-6;            %Faraday const in units [C/umol]. Note that F=Na*e, ie. charge per mole of e, elementary charge
rtf = (R*T)/F;         %in units [J/C]=V. Note that rtf is equivalent to (k*T)/e where is elementary charge and k is Boltzmann const. in [J/K]

%store constants
myC = zeros(2,1);  
myC(1) = F; myC(2) = rtf; 


%% Parameters- general 
nMit = 882;         % %number of mitochondria per cell, unitless,BioNumbersID:109409
mitA = 0.193*nMit;        %surface area of the mitochondria in units [um2], BioNumbersID:109405
mitV = 0.29*nMit;        %volume of the mitochondria in units [um3], BioNumbersID:109399
cellA = 1600;        %surface area of the plasma membrane in units [um2], BioNumbersID:103718
cellV = 3700;        %volume of the cell in units [um3], BioNumbersID:105879

%membrane potential
mPot_plsm = -50/1000;          %membrane potential over plasma membrane in [V], BioNumbersID:114263. 
mPot_mit = -20/1000;         %membrane potential over mitochondria membrane in [V], BioNumbersID:101102.

%store general parameters
p = zeros(6,1);  
p(1)=mitA; p(2)=mitV; p(3)=cellA; 
p(4)=cellV; 
p(5)= mPot_plsm; 
p(6)=mPot_mit;

%% Parameters - substrate(ion) specific

%An arbitrary Cation
pcC = 1e-1;          %permeability to cell membrane in [um/s]. BioNumbersID:112546
pcCm = 1e-1;         %permeability to mitochondria membrane in [um/s]. Used same value as pcK
zC = 1;               %charge valency, unitless
C_out = 150*1e-15;            %concentration of ion outside the cell in units [uM] converted to uM/um3
pcC_bindOn = 0;%1e-02;        %on binding in cell cytosol in [1/s].
pcC_bindOff = 0;%1e-03;       %off binding in cell cytosol in [1/s].
mt_bindOn =  0;%1e-02;              %on binding in cell mitochondria in [1/s].
mt_bindOff = 0;%1e-02;              %off binding in cell mitochondria in [1/s].
delta = 0%1e-5;                   % determines the effect of mitochondrial binding on mitochondrial membrane potential
N = 4;                         % Hill coefficent 
K = 126.85;                   % dissoication constant [uM]

%Set parameters per substrate
pS = zeros(11,1);        %a matrix of size nParameters (12) by nSpecies (1)

pS(1,1) = pcC; 
pS(2,1) = pcCm; 
pS(3,1) = zC; 
pS(4,1) = C_out; 
pS(5,1) = pcC_bindOn;
pS(6,1) = pcC_bindOff;
pS(7,1) = mt_bindOn; 
pS(8,1) = mt_bindOff; 
pS(9,1) = delta;
pS(10,1) = N;
pS(11,1) = K; 

 




