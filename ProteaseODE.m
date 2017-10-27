function [xdot, outpara] = ProteaseODE(t, y, P1)
Amino = y(1); %Amino acids
Ribo = y(2);  %Ribosomes
TC = y(3); %TC
Protein = y(4);  %Proteins
mRNA = P1.mRNAPool(1);  %mRNA Metabolic Proteins

%Values, Switch, Rates
k_codon_slow =P1.k_cs * Amino./(Amino+P1.km);
k_codon_fast = P1.k_cf * Amino./(Amino+P1.km);
Ribo_RBS = P1.RBS .* Ribo;

TC_max = P1.x_cs + (P1.L - P1.x_cs) .* (k_codon_slow./k_codon_fast);
TC_minQ = P1.x_cs .* Ribo_RBS ./k_codon_fast;
Q = (Ribo_RBS - k_codon_slow) > P1.epsi;
TC_below_max = (TC - TC_max) < -P1.epsi;
TC_above_minQ = (TC - TC_minQ) > P1.epsi;

k_init = Q.*(TC_below_max.*Ribo_RBS+(1-TC_below_max).*k_codon_slow)+(1-Q).*Ribo_RBS;

Ribo_dist = k_codon_fast./ k_init;
TC_min = P1.L ./ Ribo_dist;
TC_above_min = (TC - TC_min) > P1.epsi;
TC_equals_min = abs(TC - TC_min) < P1.epsi;

k_release = Q.*TC_above_minQ.*k_codon_slow + (1-Q).*(TC_above_min.*k_codon_slow +(1-TC_above_min).*TC_equals_min.*Ribo_RBS);

k_Acat = P1.Acat;   
k_Adep = Q.*k_codon_slow.*P1.L+(1-Q).*Ribo_RBS.*P1.L;

%k_ProteaseDeg =(Protein*P1.tag*Protease)./ (P1.Kdeg + Protein) ; %Protease degradation rate Michaelis Menten


%-------------------------------------------------------------------------------------------------------------------------------------------------
dAmino = k_Acat - k_Adep;% - dilution*Amino; %k_ProDeg
dRibo =  0;% k_ProdTotal - Total_TC  - dilution*Ribo; %k_Prod_Ribo
dTC = k_init - k_release;% - dilution*TC_Meta;
dProtein = k_release(1);% - dilution*Meta;
dProtease = k_release(2);
%-------------------------------------------------------------------------------------------------------------------------------------------------

xdot = [dAmino; %x(1)
    dRibo; %x(2)
    dTC; %x(3)
    dProtein; %x(4)
    dProtease; %x(5)
    ];
outpara = [k_init;
    k_release;
    Ribo_RBS;
    ];

end