ribo_density_max = 10;
max_density = ribo_density_max/2;
alpha = 1;
k_slow = 2;
n=1;
a=20;
ribo_density = 0:0.2:200;
for i=1:length(ribo_density)
  k_Prod(i)=k_slow.* ribo_density(i)^n ./(max_density^n + ribo_density(i)^n);
  %k_Prod(i) = (alpha.*ribo_density(i)-k_slow).*(max_density.^a ./(max_density.^a + ribo_density(i).^a))+k_slow;  
  k_Prod2(i) = k_Prod(i)./ribo_density(i);
end
figure
%plot(ribo_density,k_Init), hold on
plot(ribo_density,k_Prod),hold on
plot(ribo_density,k_Prod2)
plot(ribo_density,(k_slow)*ones(size(ribo_density)), 'LineStyle','--','LineWidth', 2)
ylabel ('Rate');
xlabel ('Ribosome Density');
legend ('Protein Production Rate', 'Prot Prod Rate/Ribo Dens','Rate slow codon');
title ('Input Graph: Protein Production Rate Sigmoid function');


ribo_density2 = 0:0.1:30;
figure
for a=2:1:20;
for i=1:length(ribo_density2)
 k_Init(i)=alpha.* ribo_density_max.^a ./(ribo_density_max.^a + ribo_density2(i).^a);
end

plot(ribo_density2,k_Init), hold on

end

legend ('a =2','a=3','a=4','a=5','a=6','a=7','a=8','a=9','a=10','a=11','a=12','a=13','a=14','a=15','a=16','a=17','a=18','a=19','a=20');
ylabel ('Rate');
xlabel ('Ribosome Density');
title ('Input Graph: Initiation Rate Sigmoid function with varied exponents');