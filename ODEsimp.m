function xdot = ODEsimp(t, y, P1)

Amino = y(1);              %Amino acids
TC_Ribo = y(2); %TC Ribo.
Ribo = y(3);  %Ribosomes
TC_Meta = y(4); %TC Meta.
Meta = y(5);  %Metabolic Proteins
TC_Syn = y(6); %TC Syn.
Syn = y(7);   %Synthetic
TC_Pro = y(8); %TC Pro.
Pro = y(9);	%Protease

mRNA_Ribo = P1.mRNAPool(1);  %mRNA Ribosomes
mRNA_Meta = P1.mRNAPool(2);  %mRNA Metabolic Proteins
mRNA_Syn = P1.mRNAPool(3);   %mRNA Synthetic
mRNA_Pro = P1.mRNAPool(4);	%mRNA Protease


%Rates
tps_cf = 1./(P1.k_cf.*Amino);
k_init = P1.RBS.*Ribo; % Initiation per mRNA
tps_avgmRNA = P1.L.*tps_cf;
k_elon = 1./tps_avgmRNA;
Total_TC = sum(P1.RBS.*P1.mRNAPool*Ribo);  %Total rate of ribosome binding

k_Prod_Meta = k_elon(1)*TC_Meta;
k_Prod_Ribo = k_elon(2)*TC_Ribo;
k_Prod_Syn = k_elon(3)*TC_Syn;
k_Prod_Pro = k_elon(4)*TC_Pro;
k_ProdTotal = k_Prod_Meta + k_Prod_Pro + k_Prod_Ribo + k_Prod_Syn; %Total protein production rate

k_Acat = (Meta * P1.vs * P1.s0)./ (P1.Ks + P1.s0);   % The rate of substrate intake from the environment.
k_Adep = k_Prod_Meta + k_Prod_Pro + k_Prod_Ribo + k_Prod_Syn;

k_ProDeg =  (Pro*P1.tag*Syn)./ (P1.Kdeg + Syn) ; %Protease degradation rate Michaelis Menten
%k_ProDeg =  10;

dilution = k_Adep/P1.CellMass;
%-------------------------------------------------------------------------------------------------------------------------------------------------
Amino =  k_Acat + k_ProDeg - k_Adep - dilution*Amino;
TC_Ribo = k_init(1)*mRNA_Ribo - k_Prod_Ribo - dilution*TC_Ribo;
Ribo = k_Prod_Ribo - Total_TC + k_ProdTotal - dilution*Ribo;
TC_Meta = k_init(2)*mRNA_Meta - k_Prod_Meta - dilution*TC_Meta;
Meta = k_Prod_Meta - dilution*Meta;
TC_Syn = k_init(3)*mRNA_Syn - k_Prod_Syn - dilution*TC_Syn;
Syn = k_Prod_Syn - k_ProDeg - dilution*Syn;
TC_Pro = k_init(4)*mRNA_Pro - k_Prod_Pro - dilution*TC_Pro;
Pro = k_Prod_Pro - dilution*Pro;
%-------------------------------------------------------------------------------------------------------------------------------------------------

xdot = [Amino; %x(1)
    TC_Meta; %x(2)
    Meta; %x(3)
    TC_Ribo; %x(4)
    Ribo; %x(5)
    TC_Syn; %x(6)
    Syn; %x(7)
    TC_Pro; %x(8)
    Pro; %x(9)
    ];


end
