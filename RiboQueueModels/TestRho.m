tic

clear
tps =  5;


P1 = struct('s0', 1000, 'mRNAPool', 1,...
	'x_cs',4, 'k_cs', 0.5, 'k_cf', 3,...
    'RBS',1,'L', 10, 'n', 20, 'a', 10, ...  
     'CellMass', 100);

%Initial values 
%y0 = ones(1,(1+2*P1.ngenes));
y0 = [2 5 1 1];

%ODE solver
[T1, x] = ode23(@(t,y) ODErho(t, y, P1), [0 tps], y0);

%Parameter Extract

for i = 1: length(T1)
[x2, paraout] = ODErho(T1(i),x(i,:),P1);
para_init(i,1) = paraout(1);
para_prod(i,1) = paraout(2);
para_Queue(i,1) = paraout(3);
para_density(i,1) = paraout(4);
para_RBSr(i,1) = paraout(5);
end

%Plots
figure
plot (T1, x(:, 1), 'k'), hold on
plot (T1, x(:, 2), 'g')
plot (T1, x(:, 3), 'r')
plot (T1, para_prod, 'c')

xlabel ('Time');
ylabel ('[]');
legend ('Amino Acids','Ribosomes', 'TR','Protein Production rate');

figure
plot (T1, para_prod, 'k'), hold on
plot (T1, para_init, 'g')
xlabel ('Time');
ylabel ('Rates');
legend ('Protein Production rate', 'Initiation Rate');


clear 
P2 = struct('tps', 20, 'var', 20);
RBS_Sweep(P2);
Time_Sweep(P2);



%%
clear
tps =  5;


P1 = struct('s0', 1000, 'mRNAPool', 1,...
	'x_cs',4, 'k_cs', 2, 'k_cf',1,...
    'RBS',1,'L', 10, 'n', 10, 'a', 1, ...  
     'CellMass', 100);

%Initial values 
%y0 = ones(1,(1+2*P1.ngenes));
y0 = [2 5 1 1];

%ODE solver
[T1, x] = ode23(@(t,y) ODErho(t, y, P1), [0 tps], y0);

%Parameter Extract

for i = 1: length(T1)
[x2, paraout] = ODEconti(T1(i),x(i,:),P1);
para_init(i,1) = paraout(1);
para_prod(i,1) = paraout(2);
para_Queue(i,1) = paraout(3);
para_density(i,1) = paraout(4);
para_RBSr(i,1) = paraout(5);
end

%Plots
figure
plot (T1, x(:, 1), 'k'), hold on
plot (T1, x(:, 2), 'g')
plot (T1, x(:, 3), 'r')
plot (T1, para_prod, 'c')
xlabel ('Time');
ylabel ('[]');
legend ('Amino Acids','Ribosomes', 'TR','Protein Production rate');

figure
plot (T1, para_prod, 'k'), hold on
plot (T1, para_init, 'g')
xlabel ('Time');
ylabel ('Rates');
legend ('Protein Production rate', 'Initiation Rate');



toc