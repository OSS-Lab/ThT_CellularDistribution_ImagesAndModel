
%driver script for the ion flows model. 
%This model involves two ions, ThT and TMRM, where flux across plasma and 
%mitochondria membrane is simulated. For ThT, binding and photosensitisation in both 
%cytoplasm and mitochondrial conpartments is also simulated where the bound, 
%photosensitised ThT in the mitochondrial can effect the mitochondrial MP. 
%The flux accounts only for passive movement under the influence of MP and conc. 

clear 
clc

%% INITIALISE

%call script that creates all constants and relevant parameters
twoIon_wPhoto_wVmeffects_constants

%Set initial ion condition. NOTE: steady state is high outside and zero in plasma. 
y0 = zeros(6,1);

%% RUN SIMS
%Specify the relative and absolute error tolerance:
options = odeset('RelTol',1e-6,'AbsTol',1e-6);

%set simulation time and run 1st simulation
simT = 5000;
int = 0.1;
simPoints = 0:int:simT;
[X,Y] =ode23s(@(t,y0)twoIon_wPhoto_wVmeffects_ode(t,y0,p,pS,myC),simPoints,y0,options); 

%Record result from previous sim as new initial conditions and introduce perturbation 
y0=Y(end,:);

mP = p(6,1); % save mitochondrial membrane potential for plotting

% Perturbation
pS(8,1) = 1e-2;   
pS(6,1) = 1e-2;

%set simulation time and run 2nd simulation
simT2 = 1080;
simPoints2 = 0:int:simT2;
 
%continue simualation
[X1,Y1] =ode23s(@(t,y0)twoIon_wPhoto_wVmeffects_ode(t,y0,p,pS,myC),simPoints2,y0,options); 

%% PLOT RESULTS

%Combine simulations
X2=[X;X1+X(end)];
Y2=[Y;Y1];

sec = 119/int; 
last = simT/int; 
start = last - sec; 
Yb = Y(start:last,:);
Yshort = [Yb;Y1];
Xb = X(start:last);
Xshort=[Xb;X1+Xb(end)];

% Normalsie 
first = Yshort(1,:); 
Yn = Yshort./ first; 

mP2 = p(6,1); 
e = length(X);
mphoto_before = Y2(1:e,4);
mphoto_after = Y2(e+1:end,4);
mtVm = (mphoto_before/mitV) * delta + mP;
mtVm_a =(mphoto_after/mitV) * delta + mP2;
mtVmf = [mtVm;mtVm_a];

% Calculate mitochondrial membrane potential 
e = length(Yshort);
mphoto_before = Yshort(1:e,4);
mphoto_after = Yshort(e+1:end,4);
mtVm = (mphoto_before/mitV) * delta + mP;
mtVm_a =(mphoto_after/mitV) * delta + mP2;
mtVmf = [mtVm;mtVm_a];

ThTp = log(Yn(:,3)) + log(Yn(:,4));

% potential normalised 
f = mtVmf(1); 
mtVmfnorm = mtVmf/f;

%plot
figure(1)
colororder({'k','k'})
yyaxis left  
line(Xshort,log(Yn(:,4)),'linewidth',2,'LineStyle','-.','color','b')
line(Xshort,log(Yn(:,6)),'linewidth',2,'color','r') 

xticks(Xshort(1):119:Xshort(end));
xlim([Xshort(1) Xshort(end)]);
xticklabels({'0','2','4','6','8','10','12','14','16','18','20'})
xlabel('Time (mins)');
ylabel('Log Normalised Amount');

yyaxis right 

plot(Xshort,mtVmfnorm,'linewidth',2,'color','k')
xlabel('Time (mins)');
legend({'ThT_p_h_o_t_o';'TMRMm';'ΔΨm'});
ylabel('Normalised Membrane Potential');
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',18); 

