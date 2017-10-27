tic

clear
tps =  50;


P1 = struct('s0', 1000, 'mRNAPool', 1,...
	'x_cs',25, 'k_cs', 5, 'k_cf', 25,...
    'RBS',5,'L', 50, 'epsi', 10.^-1);

%Initial values 
y0 = [2 1 0 0];

%ODE solver
[T1, x] = ode23(@(t,y) ProteaseODE(t, y, P1), [0 tps], y0);

%Parameter Extract

for i = 1: length(T1)
[x2, paraout] = ProteaseODE(T1(i),x(i,:),P1);
para_init(i,1) = paraout(1);
para_release(i,1) = paraout(2);
para_RBSr(i,1) = paraout(3);
end
 
%Plots
figure
plot (T1, para_init), hold on
plot (T1, para_release)
xlabel ('Time');
ylabel ('Rates');
legend ('Initiation Rate','Release Rate');

figure
plot (T1, x(:, 3)), hold on
plot (T1, x(:, 1))
%plot (T1, x(:, 4))
xlabel ('Time');
ylabel ('[]');
legend ('TC','Amino Acids'); %, 'Protein'

clear 
P2 = struct('tps', 100, 'var', 20);
Transient_Sweep(P2);

toc