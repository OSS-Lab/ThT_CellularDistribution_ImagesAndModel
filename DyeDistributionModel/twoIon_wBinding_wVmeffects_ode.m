
%function describing the ODEs for the two ion flow model
%input: general and substrate (ion) specific parameters and initial
%ion amounts in [moles]
%output: dydt for state variables

function dYdt = twoIon_wBinding_wVmeffects_ode(T, Y, p, pS, myC)

%% Retrieve constants and parameters
F = myC(1);
rtf = myC(2);

% Parameters- general
mitA = p(1);
mitV = p(2);
cellA = p(3);
cellV = p(4);
mPot_plsm = p(5);
mP = p(6);

%Parameters - ion specific
pcC = pS(1,1);
pcCm = pS(2,1);
zC = pS(3,1);
ThT_out = pS(4,1);
tmrm_out = pS(5,1);
Kon_c = pS(6,1);
Koff_c = pS(7,1);
Kon_m = pS(8,1);
Koff_m = pS(9,1);
delta = pS(10,1);
N = pS(11,1);
K = pS(12,1);
K2 = pS(13,1);

%% Retrieve current concentrations in the plasma
ThT_c = Y(1)/cellV;        %cytoplasm ThT conc in [umoles/um3]
ThT_m = Y(2)/mitV;         %mitochondrial ThT conc in [umoles/um3]
ThT_cbound = Y(3)/cellV;   %cytoplasm bound ThT conc in [umoles/um3]
ThT_mbound = Y(4)/mitV;    %mitochondrial bound ThT conc in [umoles/um3]
tmrm_c = Y(5)/cellV;       %cytoplasm tmrm conc in [umoles/um3]
tmrm_m = Y(6)/mitV;        %mitochondrial tmrm conc in [umoles/um3]

%% Calculate mitochondrial membrane potential

mPot_mit = delta * ThT_mbound + mP;

%% Calculate passive current flows in [moles/sec]
% implements Goldman current equation, with P (as devised by
% Hodgkin-Katz
%J = P*F*(1/rtf)*mPot*((Cout*exp(-z*(1/rtf)*mPot)-Cin)/(exp(-z*(1/rtf)*mPot)-1))

%passive outward current flow from cytosol to out for C in [moles/sec]
alpha = pcC*F*(1/rtf)*mPot_plsm; %has units [um.C/sec.umol]
beta = zC*(1/rtf)*mPot_plsm; %unitless
J_C_plsm = cellA*alpha*(ThT_out*exp(-1*beta)-ThT_c)/(exp(-1*beta)-1); %has units [C/sec]
J_C_plsm_Moles = J_C_plsm/F; %divide by F to convert from [C/sec] to [umoles/sec]

%passive outward current flow from cytosol to out for tmrm in [moles/sec]
J_tmrm_plsm = cellA*alpha*(tmrm_out*exp(-1*beta)-tmrm_c)/(exp(-1*beta)-1); %has units [C/sec]
J_tmrm_plsm_Moles = J_tmrm_plsm/F; %divide by F to convert from [C/sec] to [umoles/sec]

%passive outward current flow from mitochondria to cytosol for C in [moles/sec]
alphaM = pcCm*F*(1/rtf)*mPot_mit; %has units [um.C/sec.umol]
betaM = zC*(1/rtf)*mPot_mit; %unitless
J_C_mit = mitA*alphaM*(ThT_c*exp(-1*betaM)-ThT_m)/(exp(-1*betaM)-1); %has units [C/sec]
J_C_mit_Moles = J_C_mit/F; %divide by F to convert from [C/sec] to [umoles/sec]

%passive outward current flow from mitochondria to cytosol for tmrm in [moles/sec]
J_tmrm_mit = mitA*alphaM*(tmrm_c*exp(-1*betaM)-tmrm_m)/(exp(-1*betaM)-1); %has units [C/sec]
J_tmrm_mit_Moles = J_tmrm_mit/F; %divide by F to convert from [C/sec] to [umoles/sec]

K = K*mitV; %K is given uM. multiply by volume to get to umoles
K2 = K2*cellV; %K2 is given uM. multiply by volume to get to umoles

%% Define ODE for ionic movement in units [moles/sec]
%passive + active movement, outwards. i.e. we are tracking inside moles.
dYdt(size(Y,1),1) = 0;

dYdt(1) = -J_C_plsm_Moles+J_C_mit_Moles-Kon_c*(Y(1)^N/(K2^N+(Y(1)^N)))+Koff_c*Y(3); % ThT_c
dYdt(2) = -J_C_mit_Moles-Kon_m*(Y(2)^N/(K^N+(Y(2)^N)))+Koff_m*Y(4); % ThT_m
dYdt(3) = Kon_c*(Y(1)^N/(K2^N+(Y(1)^N)))-Koff_c*Y(3); % ThT_cbound
dYdt(4) = Kon_m*(Y(2)^N/(K^N+(Y(2)^N))) - Koff_m*Y(4); % ThT_mbound
dYdt(5) = -J_tmrm_plsm_Moles+J_tmrm_mit_Moles; % TMRM_c
dYdt(6) = -J_tmrm_mit_Moles; % TMRM_m


end
