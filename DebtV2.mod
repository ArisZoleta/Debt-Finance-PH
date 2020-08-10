//****&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&******
//Is Debt Finance Deficit Better than Tax Revenue (July 2020)
//Aris Zoleta DLSU Graduate School of Economic
//****&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*****

//********************
//ENDOGENOUS VARIABLES
//********************

var A C B G L W KP KG MC RM RB RP IVP  TR Y Yv Pn Pv P Y_bar MU thau phie psi psi_a ;

//**********************
//EXOGEONEOUS VARIABLES
//*********************

varexo e_a e_i e_g psi_i psi1_g;


//******************
//PARAMETERS
//******************

parameters  betha gamha hab kapha deltha alphha alphhag zeta teta tc tk tw rho_g rho_a rho_rm 
xi1_y xi1_thau xi1_b  phi_i phi_y phi_phie phie_bar;

//Calibration


betha = 0.879;		 //discount factor 
gamha = 0.5678;		//houshold relative risk avertion
hab  = 0.800;		// habit parameter 
kapha = 0.5678;		//Frisch Labor supply elasticity
deltha = 0.3456; 	//Capital depreciation rate 
alphha = 0.60;		//Share of Capital 
alphhag = 0.50;		//Share of Govt Capital
zeta = 0.56789;		//Elasticity of subs Intermediate goods 
teta = 0.678; 		//Calvo parameter
tc = 0.5609;		//Consumption tax
tk = 0.00001;		//Capital Earning tax
tw =0.00001;		//Wage earning tax 
rho_g = 0.004567;	//Govt persistence
rho_a = 0.50;		//Technology persistence 
rho_rm = 0.4589;	//Monetary policy persistence 
xi1_y = 0.6786; 	// Fiscal rule1 Output coefficient
xi1_thau = 0.890;	//Fiscal rule1  Deficit coefficient
xi1_b = 0.9000;		//Fiscal rule1  Debt coefficient
phi_i = 0.6789;		// Taylor Rule Interest rate 
phi_y = 0.4567;		//  Taylor Rule output
phi_phie = 0.6789;	// Taylor Rule Inflation
phie_bar = 0.2500;	// Target Inflation rate 

//Steady State

Ass =1;
Gss =1;
Pss = 1;
Rss = Pss*((1/betha)-(1-deltha));
MUss = (Css^(-gamha))*((1-hab)^(-gamha))*(1-hab*betha);
MCss = ((zeta-1)/zeta)*(1-betha*teta)*Pss;
Wss =  (1-alphha)*(MCss^(1/(1-alphha)))*((alphha/Rss)^(alphha/(1-alphha)));
Yss = ((Rss/(Rss-deltha*alphha))^(gamha/(gamha+kapha)))*((Wss/Pss)*(Wss/((1-alphha)*MCss))^kapha)^(1/(gamha+kapha));
Kss = alphha*MCss*(Yss/Rss);
Iss = deltha*Kss;
CRss = Yss-Iss;                  //(1/(Yss)^(kapha/gamha))*((Wss/Pss)*(Wss/(1-alphha)*MCss)^kapha)^(1/gamha);
Lss = (1-alphha)*MCss*(Yss/Wss);


//****************
//MODEL EQUATION 
//******************

//********************************
//HOUSEHOLD First Order Condition
//*********************************


//**Model Begins Here***//

model;

// Eq(1) : Household MRS

MU= (((C-hab*C(-1))^(-gamha))-hab*betha*((C(+1)-hab*C)^(-gamha)))/(1+tc);

//Eq(2) : Equation (6) and (7) Labor Suppy

W = (L^kapha)/(1+tw)*MU;

//Eq (3) : Equation (3) Law of Motion of Capital

KP = (1-deltha)*KP(-1) + IVP;

//Eq(4) : Equation (8) SDF

(1+RB) = betha*(MU(+1)/MU)*phie;

//Eq(5) : Euler Equation

1 = betha*(MU(+1)/MU)*(1-deltha + (1+tk)*RP/P);

//***************
//FIRM
//***************


//Eq (6) : Intermediate Cobb-Douglas production function 

Yv =A*(KP(-1))^alphha*(L)^1-alphha*(KG)^alphhag;

//Eq (7) Demand for Intermediate Goods 

Y = (Pv/P)^zeta*Yv; 

//Eq (8) : Firms demand for labor 

W*L=(1-alphha)*Y*MC;

//Eq (9) : Demand for Private Capital

RP*KP = alphha*Y*MC; 

//Eq (10) : Marginal Cost 

MC = (1/A*KG^alphhag)*(W/1-alphha)^1-alphha*(RP/alphha)^alphha;

//Eq (11) Optimal Price 

Pn = (zeta/(zeta-1))*(1/1-betha*teta)*MC;

//Eq (12) General Price 

P^(1-zeta) = (1-teta) * Pn^(1-zeta) + teta* P(-1)^(1-zeta) ;

//Eq (13) Price of Intermediate Goods 

Pv^(-zeta) = (1-teta)*Pn^(-zeta) + teta*P(-1)^(-zeta); 

//Eq (14) : Gross inflation

phie= P/P(+1); 

//************
//GOVERNMENT
//*************

//Eq (15) Government Budget Constraint 

B(+1)/P(+1) = (RB/phie)*B/P + thau + psi;

//Eq (16) Govt Primary Deficit 

thau = TR - G;

//Eq (17) 

TR = (1+tc)*C +(1+tk)*RP*KP +(1+tw)*W*L; 
 
//Eq (18) Govt flow of Capital 

KG = (1-deltha)*KG(-1) + G;


//*******************
//FISCAL RULES
//*******************

//Eq (20a) Rule 1 : The gov respond on output,deficit to output and debt to GDP

G = xi1_y*Y + xi1_thau*(thau/Y) + xi1_b*(B/Y) + psi1_g;

//Eq (20b) Rule 2 : The gov now repond on deficit to output and Debt to GDP relative from its target 

//G = xi2_y*Y +  xi2_thau*((thau/Y)/(thau_bar))  + xi1_b*((B/Y)/(B_Y_bar))  + psi2_g;

//*********************
//MONETARY AUTHORITY
//*********************

//Eq(21) Taylor Rule 

RM = phi_i*RM(-1) + phi_y*(Y_bar-Y) + phi_phie*(phie_bar-phie) + psi_i;

//***************
//MARKET CLEARING  
//***************


//Eq(22) Equilibrium

Y = C+ IVP+ G;

// Eq (23) Steady State output 

Y_bar = log( Y/Y(+1));

//*******************
//SHOCKS
//*****************

//Eq (24) Productivity shock

A  = rho_a*A(-1)+ e_a;

//Eq (25)

psi_i = rho_rm*psi_i(-1) + e_i;   
 
//Eq (26)     

psi1_g = rho_g*psi1_g(-1) + e_g;

end;

//***Model Ends Here***//


//*******************
//INITIAL VALUE
//******************

initval;


A = Ass;			//technology
C = Css;		//household consumption 
B = 0.90;		//government bond
G = Gss;		//Government expedinture 
L = Lss;	 	// labour hours
W = Wss;		// wages 
KP = Kss;		// private capital
KG = 0.45;		//govt capital
MC = MCss;		// marginal cost 
RM = 0.45;		// monetary policy rate 
RB = 0.30;		// bond gross interest rate 
RP = Rss;		//rental rate on capital 
IVP = Iss;		//Private Investment 
Y = Yss ;		//market clearing		//Government Investment
Yv = 0.76;		//is intermediate firm production function
Pn  = 0.6786;		//Optimal Price 
Pv  = 0.453;		//Intermediate good price
P  = Pss;		//price level/gross price 
TR = 3.2;		// Govt total revenue/ Tax collection 
MU = MUss;		// household MRS 
thau = 0.567;		//primary deficit 
phie = 0.785;		//gross inflation
psi = 4.50000;      // Govt Transfer
Y_bar = 1.56600;    //Steady State Output
psi_a = 1.008;      //productivity shock
psi_i = 1.1234;       //Monetary policy shock 
psi1_g = 0.90777;    //Government shock


end;

//Steady 

steady;


// Blanchard-Kahn conditions
check;


// Perturbation analysis
shocks;
var e_a; 
stderr 0.10;
var e_i; 
stderr 0.10;
end;


//Simulation 

simul(periods = 200);





