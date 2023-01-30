
%driver script for the ion flows model. 
%This model involves two ions, ThT and TMRM, where flux across plasma and 
%mitochondria membrane is simulated. For ThT, binding in both 
%cytoplasm and mitochondrial conpartments is also simulated where the mitochondrial bound compartment can
%effect the mitochondrial MP. The flux accounts only for passive 
%movement under the influence of MP and conc. 

%this script is used to plot TMRM steady state in the mitochondria for
%different parameter values of Kon and ThTout in order to simulate
%concentration and blue light effects on mitochondrial membrane potential. 

clear 
clc

%%
  
%call script that creates all constants and relevant parameters
twoIon_wBinding_wVmeffects_constants

n = 25; % number of simulations  
Kon = logspace(-3, -1,n); % parameter values for both Kon_c and Kon_m 
ThT_out = linspace(150e-15,1e-15,n); % parameter values for the outside concentration of ThT

result = zeros(n,n); 
count = 1;

for i = 1:length(Kon)
    
    pS(8,1) = Kon(i);
    pS(6,1) = Kon(i);
    
    for j = 1:length(ThT_out)
        
        pS(4,1) = ThT_out(j);
        
        %% RUN SIMS
        %Set initial ion condition.
        y0 = zeros(6,1);
        
        %Specify the relative and absolute error tolerance:
        options = odeset('RelTol',1e-6,'AbsTol',1e-6);
        
        %set simulation time and run 1st simulation
        simT = 2500;
        int = 10;
        simPoints = 0:int:simT;
        [X,Y] = ode23s(@(t,y0)twoIon_wBinding_wVmeffects_ode(t,y0,p,pS,myC),simPoints,y0,options);
        
        % Save TMRM_m end point
        result(j,i) = Y(end,2);
        count = count+1;
    end
end
 
% Normalise results 
normresult = zeros(n,n);

for k = 1:n
    
  norm  = result(k,:)/result(k,1); 
  normresult(k,:) = norm;
    
end 

%% Plot results 

figure (1) 
imagesc(normresult);
colorbar 

xticks(1:5:n) 
xticklabels(strsplit(num2str(Kon)))
yticks(1:5:n)
yticklabels(strsplit(num2str(ThT_out)))

ylabel('ThTout')
xlabel('Kon_c & Kon_m') 
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',20); 
