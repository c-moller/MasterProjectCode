function slowCodonSweep(P2)

tps=P2.tps;
endSweep = P2.var;
vardeci = 0.2;
for var = 1:endSweep

 P1 = struct('s0', 1000, 'mRNAPool', 1,...
	'x_cs',6, 'k_cs', vardeci, 'k_cf', 2,...
	'ribo_density_max',2,...
    'RBS',2,'L', 10, 'n', 20, 'a', 1, ...  
     'CellMass', 100);

%Initial values 
%y0 = ones(1,(1+2*P1.ngenes));
y0 = [2 4 1 1];

%ODE solver
[T1, x] = ode23(@(t,y) ODErho(t, y, P1), [0 tps], y0);
%[T1, x] = ode23(@(t,y) ODEtestT(t, y, P1), [0 tps], y0);
%[T1, x] = ode23(@(t,y) ODEconti(t, y, P1), [0 tps], y0);
%[T1, x] = ode23(@(t,y) ODEstepMin(t, y, P1), [0 tps], y0);

%Parameter Extract

for i = 1: length(T1)
[x2, paraout] = ODErho(T1(i),x(i,:),P1);
%[x2, paraout] = ODEtestT(T1(i),x(i,:),P1);
%[x2, paraout] = ODEconti(T1(i),x(i,:),P1);
%[x2, paraout] = ODEstepMin(T1(i),x(i,:),P1);

para_init(i,var) = paraout(1);
para_prod(i,var) = paraout(2);
para_Queue(i,var) = paraout(3);
para_density(i,var) = paraout(4);
para_RBSr(i,var) = paraout(5);
para_MaxDens(i,var) = paraout(6);
end
last_init(var,:)= para_init(i,var);
last_prod(var,:)= para_prod(i,var);
last_density(var,:)= para_density(i,var);
last_RBSr(var,:)= para_RBSr(i,var);
last_MaxDens(var,:) = para_MaxDens(i,var);
vardeci= vardeci+0.2;
end

% %k slow PARA SWEEP GRAPH
val=1:var;
x1plot = last_density;
x2plot = last_prod;
figure
plot(x2plot,x1plot, '* --'), hold on
plot(x2plot, para_MaxDens(1,1)*ones(size(x2plot)), 'LineStyle','--','LineWidth', 2)
plot(para_RBSr(1,1)*ones(size(x1plot)), x1plot, 'LineStyle','--','LineWidth', 2)
title('Change of Slow Codon Translation Rate effect on Ribosome density ');
xlabel('Rate of Protein Production');
ylabel('Ribosome Density per mRNA');
legend('Ribo Dens', 'Max Ribo Dens', 'RBS*Ribo');

figure
plot(val,x2plot, '* --'), hold on
plot(val, para_RBSr(1,1)*ones(size(val)),'--', 'LineWidth', 2)
plot(last_RBSr(1,1)*ones(size(x2plot)), x2plot, '--', 'LineWidth', 2)
title('Change of Slow Codon Translation Rate effect on Protein Production Rate');
xlabel('Translation Rate of Slow Codon');
ylabel('Rate Protein Production');
legend('Rate Protein Prod ','RBS*Ribo', 'RBS*Ribo');
