clear
tps=20;
y0 = [0 10 0 0];
P1 = struct('s0', 1000, 'mRNAPool', 1,...
    'x_cs',7, 'k_cs', 1, 'k_cf', 5,...
    'RBS',2,'L', 10, 'n', 20, 'a', 10, ...
    'CellMass', 100);

%ODE solver
[T1, x] = ode23(@(t,y) ODEstepMin(t, y, P1), [0 tps], y0);

%Parameter Extract

for i = 1: length(T1)
[x2, paraout] = ODEstepMin(T1(i),x(i,:),P1);
para_init(i,1) = paraout(1);
para_prod(i,1) = paraout(2);
para_Queue(i,1) = paraout(3);
para_density(i,1) = paraout(4);
para_RBSr(i,1) = paraout(5);
end

%Plots
figure
plot (T1, para_prod, 'k'), hold on
plot (T1, para_init, 'b')
plot (T1, x(:, 3), 'r')
xlabel ('Time');
ylabel ('Rates');
legend ('Protein Production rate', 'Initiation Rate', 'TC');
title({'k_{cs}=1, k_{cf}=5, mRNA length=10, Slow Codon Loc=4';'Initial values: Ribo Pool=10, TC=0'});

clear 
P2 = struct('tps', 20, 'var', 20);
RBS_Sweep(P2);
Time_Sweep(P2);