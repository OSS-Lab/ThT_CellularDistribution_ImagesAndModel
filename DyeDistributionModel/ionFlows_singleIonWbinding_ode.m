
%function describing the ODEs for the single ion flow model
%input: general and substrate (ion) specific parameters and initial
%plasma ion amounts in [moles]
%output: dydt for state variables
%assumptions: passive flux for the ion

function dYdt = ionFlows_ode(T, Y, p, pS, myC)

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
    C_out = pS(4,1); 
    pcC_bindOn = pS(5,1);
    pcC_bindOff = pS(6,1);
    mt_bindOn = pS(7,1);
    mt_bindOff = pS(8,1);
    delta = pS(9,1); 
    N = pS(10,1); 
    K = pS(11,1); 
   
    %% Retrieve current concentrations in the plasma
    C_plsm = Y(1)/cellV;    %plasma C conc in [umoles/um3]
    C_mit = Y(2)/mitV;     %mitochondria C conc in [umoles/um3]
    C_bound = Y(3)/cellV;   %plasma bound C conc in [umoles/um3]
    M_bound = Y(4)/mitV;    %mitochondrial bound C conc in [umoles/um3]
    
    %% Calculate mitochondrial membrane potential 
 
    mPot_mit = delta * M_bound + mP; 
    
    %% Calculate passive current flows in [moles/sec]
    % implements Goldman current equation, with P (as devised by
    % Hodgkin-Katz
    %J = P*F*(1/rtf)*mPot*((Cout*exp(-z*(1/rtf)*mPot)-Cin)/(exp(-z*(1/rtf)*mPot)-1))
    
    %passive outward current flow from cytosol to out for C in [moles/sec]
    alpha = pcC*F*(1/rtf)*mPot_plsm; %has units [um.C/sec.umol]
    beta = zC*(1/rtf)*mPot_plsm; %unitless
    J_C_plsm = cellA*alpha*(C_out*exp(-1*beta)-C_plsm)/(exp(-1*beta)-1); %has units [C/sec]
    J_C_plsm_Moles = J_C_plsm/F; %divide by F to convert from [C/sec] to [umoles/sec]
    
    %passive outward current flow from mitochondria to cytosol for C in [moles/sec]
    alpha = pcCm*F*(1/rtf)*mPot_mit; %has units [um.C/sec.umol]
    betaM = zC*(1/rtf)*mPot_mit; %unitless
    J_C_mit = mitA*alpha*(C_plsm*exp(-1*betaM)-C_mit)/(exp(-1*betaM)-1); %has units [C/sec]
    J_C_mit_Moles = J_C_mit/F; %divide by F to convert from [C/sec] to [umoles/sec]
    
    K = K*mitV; %K is given uM. multiply by volume to get to umoles 
        
    %% Define ODE for ionic movement in units [moles/sec]
    %passive + active movement, outwards. i.e. we are tracking inside moles.
    dYdt(size(Y,1),1) = 0; 
    dYdt(1) = -J_C_plsm_Moles+J_C_mit_Moles-pcC_bindOn*Y(1)+pcC_bindOff*Y(3);
    dYdt(2) = -J_C_mit_Moles-mt_bindOn*(Y(2)^N/(K^N+(Y(2)^N)))+mt_bindOff*Y(4); 
   % dYdt(2) = -J_C_mit_Moles-a*mt_bindOn*Y(2)+mt_bindOff*Y(4); 
    dYdt(3) = pcC_bindOn*Y(1)-pcC_bindOff*Y(3);
    dYdt(4) = mt_bindOn*(Y(2)^N/(K^N+(Y(2)^N))) - mt_bindOff*Y(4); 
   %dYdt(4) = a*mt_bindOn*Y(2) - mt_bindOff*Y(4); 
    
end 
