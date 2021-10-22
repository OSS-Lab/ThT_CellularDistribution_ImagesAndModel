
%driver script for the ion flows model. 
%This model involves only one ion, where its flux across plasma and 
%mitochondria membrane is simulated as well as ion that becomes bound in both 
%cytoplasm and mitochondrial conpartments. Mitochondrial bound compartment can
%effect the mitochondrial MP. The flux accounts only for passive 
%movement under the influence of MP and conc. 

clear 
clc

%% INITIALISE

%call script that creates all constants and relevant parameters
ionFlows_singleIonWbinding_constants

%Set initial ion condition. NOTE: steady state is high outside and zero in plasma. 
%If you start too far away from steady state, integration does not work
%column vector of size nSpecies(1) x 2 (for cell and mitochondria)
y0 = zeros(4,1);
y0(1) = (1*C_out/(10))*cellV;   % amount of ion in cell in [moles]- starting condition is away from steady state, 1/10th in cell
y0(2) = 0.0;  %mitochondrial ion in [moles]
y0(3) = 0.0;  %cytosolic bound ion in [moles]
y0(4) = 0.0; % mitochondrial bound ion in [moles] 

%% RUN SIMS
%Specify the relative and absolute error tolerance:
options = odeset('RelTol',1e-6,'AbsTol',1e-6);

%set simulation time and run 1st simulation
simT = 500000;
simPoints = 0:10:simT;
[X,Y] =ode23s(@(t,y0)ionFlows_singleIonWbinding_ode(t,y0,p,pS,myC),simPoints,y0,options); 

mP = p(6,1); % save mitochondrial membrane potential value for plotting 

%Record result from previous sim as new initial conditions and introduce perturbation 
y0=Y(end,:);

p(6) =  -75/1000; % perturbation of base mitochondrial membrane potential

%continue simualation
[X1,Y1] =ode23s(@(t,y0)ionFlows_singleIonWbinding_ode(t,y0,p,pS,myC),simPoints,y0,options); 

%combine results of two simulations into one
X1=X1+X(end);
X2=[X;X1];
Y2=[Y;Y1];

%% PLOT RESULTS
figure(1)
subplot(3,1,1)
plot(X2,Y2(:,1),'linewidth',2,'color','black');    %Cation conc. plasma
line(X2,Y2(:,3),'linewidth',2,'LineStyle','--','color','k')
xlabel('Time (sec)');
ylabel('Amount [umoles]');
legend({'ThT Plsm'; 'ThT PlsmBound'});

subplot(3,1,2)
plot(X2,Y2(:,2),'linewidth',2,'color','red');  %Cation conc. mitochondria
        line(X2,Y2(:,4),'linewidth',2,'LineStyle','--','color','red')
xlabel('Time (sec)');
ylabel('Amount [umoles]');
legend({'ThT Mit'; 'ThT MitBound'});

subplot(3,1,3)
plot(X2,Y2(:,1)/cellV,'linewidth',2,'color','red');  %Cation conc. mitochondria
        line(X2,Y2(:,2)/mitV,'linewidth',2,'LineStyle','--','color','red')
xlabel('Time (sec)');
ylabel('Conc. [umoles/um3]');
legend({'ThT Plsm'; 'ThT Mit'});

% Calculate mitochondrial membrane potential 
mP2 = p(6,1); 
delta = pS(9,1); % save delta parameter for calculating mitochondrial membrane potential

e = length(X);
mit_bound_before = Y2(1:e,4);
mitbound_after = Y2(e+1:end,4);
mtVm = (mit_bound_before/mitV) * delta + mP;
mtVm_a =(mitbound_after/mitV) * delta + mP2;
mtVmf = [mtVm;mtVm_a];
 
figure(2) 
plot(X2,mtVmf.*1000,'linewidth',2,'color','blue')
xlabel('Time (sec)');
ylabel('Membrane Potential (V)');
legend({'mPotmit'});

set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',18); 
 
