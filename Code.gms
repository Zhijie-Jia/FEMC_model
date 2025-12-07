* Definition of sets for suffix ---------------------------------------
Set     u       SAM entry     /AFAF, CEE, PFEDN, HOU, MCHS, TCI, FIN, OTH, CAP_RUR1, CAP_RUR2, CAP_RUR3, CAP_RUR4, CAP_RUR5, CAP_URB1, CAP_URB2, CAP_URB3, CAP_URB4, CAP_URB5, LAB_RUR1, LAB_RUR2, LAB_RUR3, LAB_RUR4, LAB_RUR5, LAB_URB1, LAB_URB2, LAB_URB3, LAB_URB4, LAB_URB5, IDT, TRF, DRT, SUB, RUR1, RUR2, RUR3, RUR4, RUR5, URB1, URB2, URB3, URB4, URB5, GOV, INV, DBT, ROW/
        i(u)    goods         /AFAF, CEE, PFEDN, HOU, MCHS, TCI, FIN, OTH/
        h(u)    factor        /CAP_RUR1, CAP_RUR2, CAP_RUR3, CAP_RUR4, CAP_RUR5, CAP_URB1, CAP_URB2, CAP_URB3, CAP_URB4, CAP_URB5, LAB_RUR1, LAB_RUR2, LAB_RUR3, LAB_RUR4, LAB_RUR5, LAB_URB1, LAB_URB2, LAB_URB3, LAB_URB4, LAB_URB5/
        hcap(h) CAP           /CAP_RUR1, CAP_RUR2, CAP_RUR3, CAP_RUR4, CAP_RUR5, CAP_URB1, CAP_URB2, CAP_URB3, CAP_URB4, CAP_URB5/
        hlab(h) LAB           /LAB_RUR1, LAB_RUR2, LAB_RUR3, LAB_RUR4, LAB_RUR5, LAB_URB1, LAB_URB2, LAB_URB3, LAB_URB4, LAB_URB5/
        l(u)    household     /RUR1, RUR2, RUR3, RUR4, RUR5, URB1, URB2, URB3, URB4, URB5/;
Alias (u,v), (i,j,jj), (h,k), (l,al,bl), (hcap,hcapa), (hlab,hlaba) ;
*-------------------------------------------------------------------------------

* 2. Data processing -----------------------------------------------------------
* 2.1. Loading data from SAM (excel file) --------------------------------------
*execute "gdxxrw 2018IOT_multi_HOH_group.xlsx output=data.gdx par=SAM rng=¾ÓÃñ+Õ®ÎñSAM!A1:AU47 cdim=1 rdim=1 par=additional rng=addition!B1:K2 cdim=1 rdim=0 par=epsilonLD rng=addition!B3:K4 cdim=1 rdim=0 ";
Parameter        SAM(u,v)              "Social accounting matrix"
                 additional(l)         "Outstanding debt"
                 epsilonLD(hlab)       "Labor supply wage elasticity for labor type hlab (labor supply curve)"

;
$gdxin   data.gdx
$loaddc  SAM
$loaddc  additional
$loaddc  epsilonLD
$gdxin
;

* Loading the initial values ------------------------------------------
Parameter
                 F0(h,j)              "Benchmark factor demand of factor h by activity j (from SAM)"
                 LAB0(j)              "Benchmark total labor input by activity j"
                 CAP0(j)              "Benchmark total capital input by activity j"
                 VA0(j)               "Benchmark value added by activity j"
                 X0(i,j)              "Benchmark intermediate input of good i by activity j"
                 TX0(j)               "Benchmark total intermediate demand of activity j"
                 Z0(j)                "Benchmark gross output of activity j"
                 Xp0(i,l)             "Benchmark household l consumption of good i"
                 Xg0(i)               "Benchmark government consumption of good i"
                 Xv0(i)               "Benchmark investment demand for good i"
                 E0(i)                "Benchmark exports of good i"
                 M0(i)                "Benchmark imports of good i"
                 Q0(i)                "Benchmark Armington composite quantity of good i"
                 D0(i)                "Benchmark domestic supply of good i"
                 Loan0(l)             "Benchmark new borrowing (loan inflow) of household l"
                 Loang0               "Benchmark new borrowing (loan inflow) of government"
                 refund0(l)           "Benchmark loan repayment (principal+interest) by household l"
                 refundg0             "Benchmark loan repayment (principal+interest) by government"
                 Outstanding_debt(l)  "Outstanding stock of debt of household l at benchmark"
                 Outstanding_debt_g   "Outstanding stock of government debt at benchmark"

                 HOHsalary0(l)        "Benchmark factor income (labor+capital) of household l"
                 HOHincome0(l)        "Benchmark total income of household l (factor income + net loans)"
                 Tp0(l)               "Benchmark government transfers to household l"

                 Sp0(l)               "Benchmark private saving of household l"
                 Sg0                  "Benchmark government saving"
                 deficit              "Exogenous fiscal deficit (policy shock size)"
                 Td0(l)               "Benchmark direct income tax paid by household l"
                 Tz0(j)               "Benchmark production (indirect) tax paid by activity j"
                 Tm0(i)               "Benchmark import tariff revenue on good i"
                 GOVincome0           "Benchmark total government income"
                 GOVtax0              "Benchmark total tax revenue (direct + indirect + tariffs)"
                 GDP0                 "Benchmark (real) GDP at base prices"
                 nGDP0                "Benchmark nominal GDP at current prices"
                 PPI0(j)              "Benchmark producer price index for activity j (base =100)"
                 CPI0                 "Benchmark consumer price index (base =100)"
                 Velocity0            "Benchmark income velocity of money"
                 coeff_vel            "Coefficient linking velocity to consumption share in income"
                 M2                   "Exogenous money supply (M2) at benchmark"
                 mmultiplier          "Money multiplier from base money to M2"

                 FF0(h)               "Benchmark factor endowment/total supply of factor h"
                 Debt0                "Benchmark net increase in total debt (flow)"
                 Sf                   "Foreign saving (current account deficit, in foreign currency)"
                 pWe(i)               "World export price of good i in foreign currency"
                 pWm(i)               "World import price of good i in foreign currency"
                 tauz(i)              "Production tax rate on output of good i"
                 taum(i)              "Import tariff rate on good i"

                 Investment           "Exogenous additional government investment (policy variable)"
                 LSS(l)               "Lump-sum subsidy from government to household l"
                 coupon               "Total nominal value of consumption coupon (voucher) program"
                 share_coupon(i,l)    "Share of coupon allocated to good i and household l"
;
Td0(l)  =SAM("DRT",l);
Tz0(j)  =SAM("IDT",j);
Tm0(j)  =SAM("TRF",j);

F0(h,j) =SAM(h,j);
LAB0(j) =sum(hlab,SAM(hlab,j));
CAP0(j) =sum(hcap,SAM(hcap,j));
VA0(j)   =sum(h, SAM(h,j));
X0(i,j) =SAM(i,j);
TX0(j)  =sum(i, X0(i,j));
Z0(j)   =VA0(j) +TX0(j);
M0(i)   =SAM("ROW",i);

tauz(j) =Tz0(j)/Z0(j);
taum(j)$(M0(j))=  Tm0(j)/M0(j);

Xp0(i,l)  =SAM(i,l);
FF0(h)    =sum(j,F0(h,j));

Loan0(l)         =       SAM(l,"DBT");
Loang0           =       SAM("GOV","DBT");
refund0(l)       =       SAM("DBT",l);
refundg0         =       SAM("DBT","GOV");
Debt0            =       SAM("DBT","INV");
Outstanding_debt(l)  =       additional(l);
Outstanding_debt_g   =       90100;

GOVincome0       =  sum(l, Td0(l)) +sum(j, Tz0(j)) +sum(j, Tm0(j)) +loang0;
GOVtax0          =  sum(l, Td0(l)) +sum(j, Tz0(j)) +sum(j, Tm0(j));

HOHsalary0("RUR1")  = SAM("RUR1","LAB_RUR1") + SAM("RUR1","CAP_RUR1");
HOHsalary0("RUR2")  = SAM("RUR2","LAB_RUR2") + SAM("RUR2","CAP_RUR2");
HOHsalary0("RUR3")  = SAM("RUR3","LAB_RUR3") + SAM("RUR3","CAP_RUR3");
HOHsalary0("RUR4")  = SAM("RUR4","LAB_RUR4") + SAM("RUR4","CAP_RUR4");
HOHsalary0("RUR5")  = SAM("RUR5","LAB_RUR5") + SAM("RUR5","CAP_RUR5");
HOHsalary0("URB1")  = SAM("URB1","LAB_URB1") + SAM("URB1","CAP_URB1");
HOHsalary0("URB2")  = SAM("URB2","LAB_URB2") + SAM("URB2","CAP_URB2");
HOHsalary0("URB3")  = SAM("URB3","LAB_URB3") + SAM("URB3","CAP_URB3");
HOHsalary0("URB4")  = SAM("URB4","LAB_URB4") + SAM("URB4","CAP_URB4");
HOHsalary0("URB5")  = SAM("URB5","LAB_URB5") + SAM("URB5","CAP_URB5");

HOHincome0(l)       = HOHsalary0(l) + Loan0(l);

Tp0(l)              = SAM(l,"SUB");

Xg0(i)  =SAM(i,"GOV");
Xv0(i)  =SAM(i,"INV");
E0(i)   =SAM(i,"ROW");
Q0(i)   =sum(l,Xp0(i,l))+Xg0(i)+Xv0(i)+sum(j, X0(i,j));
D0(i)   =(1+tauz(i))*Z0(i)-E0(i);
Sp0(L)  =SAM("INV",l);
Sg0     =SAM("INV","GOV");
deficit =0;
Sf      =SAM("INV","ROW");

GDP0    =  sum(i,Xv0(i)+Xg0(i)+sum(l,Xp0(i,l))+E0(i)-M0(i));
nGDP0   =  sum(i,Xv0(i)+Xg0(i)+sum(l,Xp0(i,l))+E0(i)-M0(i));
PPI0(j) =  100;
CPI0    =  100;
*2018 M2=179.2928002 tillion RMB
M2               =  179292.8002;
*Money multiplier
mmultiplier      =  5.35;
Velocity0        =  GDP0*CPI0/M2;
coeff_vel        =  Velocity0 / sum(l,( sum(i, Xp0(i,l)* 1 ) / HOHincome0(l)  ) );

pWe(i)   =       1;
pWm(i)   =       1;

Investment=      0;
LSS(l)   =       0;
coupon   =       0;
share_coupon(i,l)= Xp0(i,l)/sum((j,al), Xp0(j,al));

* Calibration ---------------------------------------------------------
Parameter
                rhou(l)         "CES parameter for household l consumption (rho = 1 - 1/elasticity)"
                rhoz(j)         "CES parameter for top-level production (VA vs intermediates)"
                rhova(j)        "CES parameter for value-added nest (labor vs capital)"
                rholab(j)       "CES parameter for labor nest (across labor types) in activity j"
                rhocap(j)       "CES parameter for capital nest (across capital types) in activity j"
                rhoinc(l)       "CES parameter for household l income composition (labor vs loans)"
                rhoginc         "CES parameter for government income composition (tax vs borrowing)"
                sigma(i)        "Elasticity of substitution in Armington aggregation for good i"
                psi(i)          "Elasticity of transformation between exports and domestic supply"
                eta(i)          "Transformed Armington parameter: (sigma - 1)/sigma"
                phi(i)          "Transformed CET parameter: (psi + 1)/psi"
;

* rho = 1-1/elasticity
* Low substitution across broad consumption categories
rhou(l)      =  1-1/0.7;
rhoz(j)      =  1-1/1.5;
rhova(j)     =  1-1/1.5;
rholab(j)    =  1-1/0.8;
rhocap(j)    =  1-1/1.8;

sigma(i)=2;
psi(i)  =2;
eta(i)=(sigma(i)-1)/sigma(i);
phi(i)=(psi(i)+1)/psi(i);
rhoinc(l)    =  1+1/0.2;
rhoginc      =  1+1/0.5;

Parameter
                deltaz(j)         "CES share parameter for VA vs intermediates in production of j"
                deltava(j)        "CES share parameter for labor vs capital in VA of j"
                deltalab(hlab,j)  "CES share parameter for labor type hlab in labor nest of j"
                deltacap(hcap,j)  "CES share parameter for capital type hcap in capital nest of j"
                alphaz(j)         "Scale parameter of top-level production CES function for activity j"
                alphava(j)        "Scale parameter of value-added CES nest for activity j"
                alphalab(j)       "Scale parameter of labor CES nest for activity j"
                alphacap(j)       "Scale parameter of capital CES nest for activity j"
                ax(i,j)           "Intermediate input coefficient of good i per unit TX(j)"

                deltainc(l)       "CES share parameter of labor income in household l total income"
                alphainc(l)       "Scale parameter of CES household income function for l"
                deltaginc         "CES share parameter of tax in government income"
                alphaginc         "Scale parameter of CES government income function"

                alpha(i,l)        "Budget share parameter of good i in household l utility function"

                mu(i)             "Share parameter of good i in government consumption"
                lambda(i)         "Share parameter of good i in investment demand"
                deltam(i)         "CES share parameter of imports in Armington composite for i"
                deltad(i)         "CES share parameter of domestic good in Armington composite for i"
                gamma(i)          "Scale parameter of Armington (CES) aggregation for good i"
                xid(i)            "CET share parameter for domestic supply in total output of i"
                xie(i)            "CET share parameter for exports in total output of i"
                theta(i)          "Scale parameter of CET transformation function for good i"
                ssp(l)            "Average propensity to save out of household l disposable income"
                ssg               "Average propensity to save out of government tax revenue"
                taud(l)           "Direct income tax rate of household l"
                tautp(l)          "Transfer payment rate (transfers as share of GDP) for household l"
                taurf(l)          "Repayment rate on outstanding private debt of household l"
                taurfg            "Repayment rate on outstanding government debt"

                UU0(l)            "Benchmark utility level of household l (from initial consumption)"
                SW0               "Benchmark social welfare (sum of utilities) at baseline"
                VV0(l)            "Benchmark money-metric utility (equivalent consumption) of l"
                V_source0(l)      "Benchmark utility index using base prices (source decomposition)"
                V_use0(l)         "Benchmark utility index using actual prices (use decomposition)"
                EV0(l)            "Benchmark equivalent variation (0 at base) for household l"
                EV_source0(l)     "Benchmark EV decomposition from income/source side of l"
                EV_use0(l)        "Benchmark EV decomposition from price/use side of l"
;
deltaz(j)       =        VA0(j)**(1-rhoz(j))/(VA0(j)**(1-rhoz(j))+TX0(j)**(1-rhoz(j)));
deltava(j)      =        LAB0(j)**(1-rhova(j))/(LAB0(j)**(1-rhova(j))+CAP0(j)**(1-rhova(j)));
deltalab(hlab,j)=        F0(hlab,j)**(1-rholab(j))/sum(hlaba, F0(hlaba,j)**(1-rholab(j)));
deltacap(hcap,j)=        F0(hcap,j)**(1-rhocap(j))/sum(hcapa, F0(hcapa,j)**(1-rhocap(j)));

alphaz(j)       =        Z0(j)/(deltaz(j)*VA0(j)**rhoz(j)+(1-deltaz(j))*TX0(j)**rhoz(j))**(1/rhoz(j));
alphava(j)      =        VA0(j)/(deltava(j)*LAB0(j)**rhova(j)+(1-deltava(j))*CAP0(j)**rhova(j))**(1/rhova(j));
alphalab(j)     =        LAB0(j)/sum(hlaba, deltalab(hlaba,j)*F0(hlaba,j)**rholab(j) )**(1/rholab(j));
alphacap(j)     =        CAP0(j)/sum(hcapa, deltacap(hcapa,j)*F0(hcapa,j)**rhocap(j) )**(1/rhocap(j));

alpha(i,l)      =        Xp0(i,l)**(1-rhou(l))/sum(j, Xp0(j,l)**(1-rhou(l)) );

ax(i,j)        =  X0(i,j)/TX0(j);
mu(i)   =Xg0(i)/sum(j, Xg0(j));
lambda(i)=Xv0(i)/(sum(l,Sp0(l)) +Sg0 +Sf -debt0);

deltainc(l)    =   HOHsalary0(l)**(1-rhoinc(l))/(HOHsalary0(l)**(1-rhoinc(l)) + loan0(l)**(1-rhoinc(l)));
alphainc(l)    =   HOHincome0(l) /(deltainc(l)*HOHsalary0(l)**rhoinc(l)+(1-deltainc(l))*loan0(l)**rhoinc(l))**(1/rhoinc(l));

deltaginc      =   GOVtax0**(1-rhoginc)/(GOVtax0**(1-rhoginc) + loang0**(1-rhoginc));
alphaginc      =   GOVincome0 /(deltaginc*GOVtax0**rhoginc+(1-deltaginc)*loang0**rhoginc)**(1/rhoginc);

deltam(i)$(M0(i))=(1+taum(i))*M0(i)**(1-eta(i))
                  /((1+taum(i))*M0(i)**(1-eta(i)) +D0(i)**(1-eta(i)));
deltad(i)$(M0(i))=D0(i)**(1-eta(i))
                  /((1+taum(i))*M0(i)**(1-eta(i)) +D0(i)**(1-eta(i)));
gamma(i)$(M0(i)) =Q0(i)/(deltam(i)*M0(i)**eta(i)+deltad(i)*D0(i)**eta(i))
                  **(1/eta(i));

xie(i)$(E0(i))   =E0(i)**(1-phi(i))/(E0(i)**(1-phi(i))+D0(i)**(1-phi(i)));
xid(i)$(E0(i))   =D0(i)**(1-phi(i))/(E0(i)**(1-phi(i))+D0(i)**(1-phi(i)));
theta(i)$(E0(i)) =Z0(i)
                  /(xie(i)*E0(i)**phi(i)+xid(i)*D0(i)**phi(i))**(1/phi(i));

ssp(l)  = Sp0(l)/( HOHsalary0(l) + Tp0(l) );
ssg     = Sg0/(sum(l,Td0(l)) + sum(j, Tz0(j) + Tm0(j)));
taud(l) = Td0(l)/HOHsalary0(l);
tautp(l)= Tp0(l)/GDP0;
taurf(l)= refund0(l)/Outstanding_debt(l);
taurfg  = refundg0  /Outstanding_debt_g ;

UU0(l)  = sum(i, alpha(i,l)*Xp0(i,l)**rhou(l))**(1/rhou(l));
SW0     = sum(l,UU0(l));
VV0(l)         =     UU0(l);
V_source0(l)  =     UU0(l);
V_use0(l)     =     UU0(l);
EV0(l)        =     0;
EV_source0(l) =     0;
EV_use0(l)    =     0;
* ---------------------------------------------------------------------

* Defining model system -----------------------------------------------
Variable
                F(h,j)          "Demand for factor h by activity j"
                LAB(j)          "Aggregate labor composite used by activity j"
                CAP(j)          "Aggregate capital composite used by activity j"
                VA(j)           "Value added of activity j"
                X(i,j)          "Intermediate demand for good i by activity j"
                TX(j)           "Total intermediate input of activity j"
                Z(j)            "Gross output of activity j"
                Xp(i,l)         "Household l consumption of good i"
                Xg(i)           "Government consumption demand for good i"
                Xv(i)           "Investment demand for good i"
                E(i)            "Exports of good i"
                M(i)            "Imports of good i"
                Q(i)            "Armington composite demand for good i"
                D(i)            "Domestic supply of good i (used domestically)"

                loan(l)         "New borrowing (loan inflow) of household l"
                refund(l)       "Debt service payment by household l"
                loang           "New borrowing of government"
                refundg         "Debt service payment by government"
                debt            "Total net increase in debt (private + government)"

                FF(h)           "Total available factor endowment/supply of factor h"

                pF(h)           "Price (wage/rental) of factor h"
                pLAB(j)         "Price of aggregate labor composite in activity j"
                pCAP(j)         "Price of aggregate capital composite in activity j"
                pVA(j)          "Price of value added in activity j"
                pTX(j)          "Price index of intermediate input bundle in activity j"
                pZ(j)           "Producer price of output of activity j"
                pq(i)           "Price of Armington composite good i"
                pq_hoh(i,l)     "Effective consumer price of good i for household l (after coupons)"
                pe(i)           "Export price of good i in domestic currency"
                pm(i)           "Import price of good i in domestic currency"
                pd(i)           "Domestic basic price of good i"
                ploan(l)        "Price of loans to household l (1 + interest rate)"
                ploang          "Price of government loans (1 + interest rate on gov debt)"
                epsilon         "Nominal exchange rate (domestic per unit of foreign currency)"

                HOHsalary(l)    "Factor income (labor+capital) received by household l"
                HOHincome(l)    "Total income of household l (factor income + net borrowing)"
                Tp(l)           "Government transfers to household l"

                Sp(l)           "Private saving of household l"
                Sg              "Government saving"
                Td              "Total direct tax revenue (sum over households)"
                Tz(j)           "Production tax revenue from activity j"
                Tm(i)           "Import tariff revenue on good i"
                GOVincome       "Total government income (taxes + borrowing)"
                GOVtax          "Total tax revenue (direct + indirect + tariffs)"

                GDP             "Real GDP"
                nGDP            "Nominal GDP"
                GDPCHK          "GDP accounting check (should be zero if model is consistent)"
                PPI(j)          "Producer price index for activity j"
                CPI             "Consumer price index"
                Velocity        "Velocity of money"
                Walras          "Walrasian slack variable (to close excess market conditions)"

                UU              "Aggregate utility index (fictitious, for welfare)"
                SW              "Social welfare = sum of household utilities"

                VV(l)           "Money-metric utility of household l (current prices)"
                V_source(l)     "Decomposition: utility index driven by income/source changes"
                V_use(l)        "Decomposition: utility index driven by price/use changes"
                EV(l)           "Equivalent variation for household l (welfare change, %)"
                EV_source(l)    "EV component from income/source changes"
                EV_use(l)       "EV component from price/use changes"
;
Equation
                 eqZ1(j)            "Top-level CES production function: output Z(j) from VA and TX"
                 eqZ2(j)            "First-order condition: relative price between VA and TX"
                 eqZ3(j)            "Zero-profit condition for activity j (value of output = inputs)"
                 eqX(i,j)           "Intermediate input demand: X(i,j) proportional to TX(j)"
                 eqpzs(j)           "Intermediate bundle price pTX(j) as cost-weighted pq(i)"
                 eqCESva1(j)        "CES value-added nest: VA(j) from LAB(j) and CAP(j)"
                 eqCESva2(j)        "First-order condition: relative price of LAB vs CAP"
                 eqCESva3(j)        "Zero-profit condition for value-added nest"

                 eqLAB1(j)          "Labor nest: composite LAB(j) from individual labor types"
                 eqLAB2(hlab,j)     "Demand for labor type hlab in activity j (CES allocation)"
                 eqCAP1(j)          "Capital nest: composite CAP(j) from capital types"
                 eqCAP2(hcap,j)     "Demand for capital type hcap in activity j (CES allocation)"

                 eqHOHsalary1(l)    "Factor income definition for rural household 1"
                 eqHOHsalary2(l)    "Factor income definition for rural household 2"
                 eqHOHsalary3(l)    "Factor income definition for rural household 3"
                 eqHOHsalary4(l)    "Factor income definition for rural household 4"
                 eqHOHsalary5(l)    "Factor income definition for rural household 5"
                 eqHOHsalary6(l)    "Factor income definition for urban household 1"
                 eqHOHsalary7(l)    "Factor income definition for urban household 2"
                 eqHOHsalary8(l)    "Factor income definition for urban household 3"
                 eqHOHsalary9(l)    "Factor income definition for urban household 4"
                 eqHOHsalary10(l)   "Factor income definition for urban household 5"

                 eqHOHincome(l)     "CES income aggregation: HOHincome(l) from salary and loans"
                 eqloan(l)          "Optimal household borrowing demand as function of income and rate"
                 eqploan(l)         "Implied salary component in income CES (normalization)"
                 eqloang            "Government borrowing rule (exogenous deficit added)"
                 eqploang           "Normalization of government loan price (set to 1)"
                 eqrefund(l)        "Debt service payments of households tied to outstanding debt"
                 eqrefundg          "Government debt service payment tied to gov outstanding debt"
                 eqdebt             "Total net increase in debt (gov + households)"

                 eqTd(l)            "Direct income tax revenue from household l"
                 eqTz(j)            "Production tax revenue from activity j"
                 eqTm(i)            "Import tariff revenue on good i"
                 eqXg(i)            "Government demand for good i (Cobb-Douglas over goods)"
                 eqGOVincome        "Government budget: definition of total income"
                 eqGOVtax           "Definition of total tax revenue (GOVtax)"

                 eqXv(i)            "Investment demand for good i (Cobb-Douglas over goods)"
                 eqSp               "Private saving function (aggregate from households)"
                 eqSg               "Government saving as fixed share of tax revenue"
                 eqTp(l)            "Transfer payments to household l (GDP-based + LSS)"

                 eqXp(i,l)          "Household l demand for good i (CES utility maximization)"
                 eqpq_hoh(i,l)      "Effective consumer price after coupons for good i, household l"
                 eqpe(i)            "Export price in local currency (via exchange rate)"
                 eqpm(i)            "Import price in local currency (via exchange rate)"
                 eqepsilon          "Balance of payments: foreign exchange market clearing"

                 eqpqs(i)           "Armington CES aggregation: Q(i) from M(i) and D(i)"
                 eqM(i)             "Import demand function from Armington CES"
                 eqD(i)             "Domestic demand derived from Armington CES"
                 eqpqs1(i)          "Special case: no imports ¡ú Q(i) = D(i)"
                 eqpqs2(i)          "Special case: no imports ¡ú pq(i) = pd(i)"
                 eqpqs3(i)          "Special case: no imports ¡ú M(i) = 0"

                 eqpzd(i)           "CET transformation: Z(i) from E(i) and D(i)"
                 eqDs(i)            "Domestic supply from CET transformation"
                 eqE(i)             "Export supply from CET transformation"
                 eqpzd1(i)          "Special case: no exports ¡ú Z relates only to D"
                 eqpzd2(i)          "Special case: no exports ¡ú pz(i) = pd(i)"
                 eqpzd3(i)          "Special case: no exports ¡ú E(i) = 0"

                 eqpqd(i)           "Market clearing for Armington composite Q(i)"
                 eqpf(h)            "Factor market clearing condition for factor h"

                 eqGDP              "Definition of real GDP (expenditure side)"
                 eqnGDP             "Definition of nominal GDP"
                 eqGDPCHK           "GDP consistency check (income-expenditure identity)"
                 eqPPI(j)           "Producer price index definition for activity j"
                 eqCPI              "Consumer price index definition"
                 eqQuantityMoney    "Quantity theory: money demand (CPI*GDP = money*velocity)"
                 eqVelocity         "Definition of money velocity (tied to consumption share)"

                 eqUU(l)            "CES utility index for household l"
                 eqSW               "Social welfare as sum of household utilities"

                 eqVV(l)            "Money-metric utility VV(l) from consumption expenditure"
                 eqV_source(l)      "Utility index focusing on income/source change (price=1)"
                 eqV_use(l)         "Utility index focusing on price/use change (income fixed)"

                 eqEV(l)            "Equivalent variation (EV) for household l in %"
                 eqEV_source(l)     "EV contribution from income/source change in %"
                 eqEV_use(l)        "EV contribution from use/price change (EV - EV_source)"

                 eqlabordemand(hlab) "Upward-sloping labor supply: wage¨Clabor supply elasticity"
;
*[domestic production] ----
eqZ1(j)..                                 Z(j)  =e=  alphaz(j)*(deltaz(j)*VA(j)**rhoz(j)+(1-deltaz(j))*TX(j)**rhoz(j))**(1/rhoz(j));
eqZ2(j)..                        pva(j)/ptx(j)  =e=  deltaz(j)/(1-deltaz(j))*(TX(j)/VA(j))**(1-rhoz(j));
eqZ3(j)..                           pz(j)*Z(j)  =e=  pva(j)*VA(j)+ptx(j)*TX(j);

eqX(i,j)..                              X(i,j)  =e=  ax(i,j)*TX(j);
eqpzs(j)..                              ptx(j)  =e=  sum(i, ax(i,j)*pq(i));

eqCESva1(j)..                            VA(j)  =e=  alphava(j)*(deltava(j)*LAB(j)**rhova(j)+(1-deltava(j))*CAP(j)**rhova(j))**(1/rhova(j));
eqCESva2(j)..                  pLAB(j)/pCAP(j)  =e=  deltava(j)/(1-deltava(j))*(CAP(j)/LAB(j))**(1-rhova(j));
eqCESva3(j)..                     pVA(j)*VA(j)  =e=  pLAB(j)*LAB(j)+PCAP(j)*CAP(j);

eqLAB1(j)..                             LAB(j)  =e=  alphalab(j)*sum(hlab, deltalab(hlab,j)*F(hlab,j)**rholab(j) )**(1/rholab(j));
eqLAB2(hlab,j)..                     F(hlab,j)  =e= (alphalab(j)**rholab(j)*deltalab(hlab,j)*pLAB(j)/ pF(hlab) )**(1/(1-rholab(j)))*LAB(j);

eqCAP1(j)..                             CAP(j)  =e=  alphacap(j)*sum(hcap, deltacap(hcap,j)*F(hcap,j)**rhocap(j) )**(1/rhocap(j));
eqCAP2(hcap,j)..                     F(hcap,j)  =e= (alphacap(j)**rhocap(j)*deltacap(hcap,j)*pCAP(j)/ pF(hcap) )**(1/(1-rhocap(j)))*CAP(j);

eqHOHsalary1("RUR1")..                  HOHsalary("RUR1")  =e= sum(j, F("CAP_RUR1",j)*pf("CAP_RUR1") + F("LAB_RUR1",j)*pf("LAB_RUR1"));
eqHOHsalary2("RUR2")..                  HOHsalary("RUR2")  =e= sum(j, F("CAP_RUR2",j)*pf("CAP_RUR2") + F("LAB_RUR2",j)*pf("LAB_RUR2"));
eqHOHsalary3("RUR3")..                  HOHsalary("RUR3")  =e= sum(j, F("CAP_RUR3",j)*pf("CAP_RUR3") + F("LAB_RUR3",j)*pf("LAB_RUR3"));
eqHOHsalary4("RUR4")..                  HOHsalary("RUR4")  =e= sum(j, F("CAP_RUR4",j)*pf("CAP_RUR4") + F("LAB_RUR4",j)*pf("LAB_RUR4"));
eqHOHsalary5("RUR5")..                  HOHsalary("RUR5")  =e= sum(j, F("CAP_RUR5",j)*pf("CAP_RUR5") + F("LAB_RUR5",j)*pf("LAB_RUR5"));
eqHOHsalary6("URB1")..                  HOHsalary("URB1")  =e= sum(j, F("CAP_URB1",j)*pf("CAP_URB1") + F("LAB_URB1",j)*pf("LAB_URB1"));
eqHOHsalary7("URB2")..                  HOHsalary("URB2")  =e= sum(j, F("CAP_URB2",j)*pf("CAP_URB2") + F("LAB_URB2",j)*pf("LAB_URB2"));
eqHOHsalary8("URB3")..                  HOHsalary("URB3")  =e= sum(j, F("CAP_URB3",j)*pf("CAP_URB3") + F("LAB_URB3",j)*pf("LAB_URB3"));
eqHOHsalary9("URB4")..                  HOHsalary("URB4")  =e= sum(j, F("CAP_URB4",j)*pf("CAP_URB4") + F("LAB_URB4",j)*pf("LAB_URB4"));
eqHOHsalary10("URB5")..                 HOHsalary("URB5")  =e= sum(j, F("CAP_URB5",j)*pf("CAP_URB5") + F("LAB_URB5",j)*pf("LAB_URB5"));

*[government behavior] ----
eqTd(l)..                                Td(l)   =e= taud(l)*HOHsalary(l);
eqTz(j)..                                Tz(j)   =e= tauz(j)*pz(j)*Z(j);
eqTm(i)..                                Tm(i)   =e= taum(i)*pm(i)*M(i);
eqXg(i)..                                Xg(i)   =e= mu(i)*(GOVincome -Sg -sum(l,Tp(l)) -refundg -coupon - Investment)/pq(i);
eqTp(l)..                                Tp(l)   =e= tautp(l)*GDP + LSS(l);

eqGOVincome..                        GOVincome   =e= sum(l, Td(l)) +sum(j, Tz(j)) +sum(j, Tm(j)) +loang*ploang;
eqGOVtax..                           GOVtax      =e= sum(l, Td(l)) +sum(j, Tz(j)) +sum(j, Tm(j));

*[investment behavior] ----
eqXv(i)..       Xv(i)   =e= lambda(i)*(sum(l,Sp(l)) +Sg +epsilon*Sf -debt +Investment  )/pq(i);


*[Finance behavior] ----
eqrefund(l)..   refund(l) =e= taurf(l)  * Outstanding_debt(l);
eqrefundg..     refundg   =e= taurfg    * Outstanding_debt_g;

eqHOHincome(l)..    HOHincome(l)  =e=  alphainc(l)*(deltainc(l)*HOHsalary(l)**rhoinc(l)+(1-deltainc(l))*loan(l)**rhoinc(l))**(1/rhoinc(l));
eqloan(l)..              loan(l)  =e=  (alphainc(l)**rhoinc(l)*(1-deltainc(l))*1/ploan(l))**(1/(1-rhoinc(l)))*HOHincome(l);
eqploan(l)..        HOHsalary(l)  =e=  (alphainc(l)**rhoinc(l)*   deltainc(l) *1/1       )**(1/(1-rhoinc(l)))*HOHincome(l);
eqloang..                   loang  =e=  loang0 + deficit;
eqploang..                 ploang  =e=  1;

eqdebt..        debt      =e= loang*ploang -refundg + sum(l, loan(l)*ploan(l) -refund(l));

*[savings] ----------------
eqSp(l)..       Sp(l)   =e= ssp(l)*( HOHsalary(l) + Tp(l) );
eqSg..          Sg      =e= ssg * GOVtax;

*[household consumption] --
eqXp(i,l)..      Xp(i,l)   =e= (alpha(i,l)/pq_hoh(i,l))**(1/(1-rhou(l)))*(HOHincome(l) -Sp(l) -Td(l) +Tp(l) -refund(l))/sum(j, alpha(j,l)**(1/(1-rhou(l)))*pq_hoh(j,l)**(rhou(l)/(rhou(l)-1)) );
eqpq_hoh(i,l)..  pq_hoh(i,l)*Xp(i,l) =e=  pq(i)*Xp(i,l) -coupon*share_coupon(i,l);

*[international trade] ----
eqpe(i)..       pe(i)   =e= epsilon*pWe(i);
eqpm(i)..       pm(i)   =e= epsilon*pWm(i);
eqepsilon..     sum(i, pWe(i)*E(i)) +Sf
                        =e= sum(i, pWm(i)*M(i));

*[Armington function] -----
eqpqs(i)$(M0(i)<>0)..      Q(i)    =e= gamma(i)*(deltam(i)*M(i)**eta(i)+deltad(i)
                            *D(i)**eta(i))**(1/eta(i));
eqM(i)$(M0(i)<>0)..        M(i)    =e= (gamma(i)**eta(i)*deltam(i)*pq(i)
                            /((1+taum(i))*pm(i)))**(1/(1-eta(i)))*Q(i);
eqD(i)$(M0(i)<>0)..        D(i)    =e= (gamma(i)**eta(i)*deltad(i)*pq(i)/pd(i))
                            **(1/(1-eta(i)))*Q(i);
eqpqs1(i)$(M0(i)=0)..                     Q(i)  =e=  D(i);
eqpqs2(i)$(M0(i)=0)..                    pq(i)  =e=  pd(i);
eqpqs3(i)$(M0(i)=0)..                     M(i)  =e=  0;

*[transformation function] -----
eqpzd(i)$(E0(i)<>0)..      Z(i)    =e= theta(i)*(xie(i)*E(i)**phi(i)+xid(i)
                            *D(i)**phi(i))**(1/phi(i));
eqE(i)$(E0(i)<>0)..        E(i)    =e= (theta(i)**phi(i)*xie(i)*(1+tauz(i))*pz(i)
                            /pe(i))**(1/(1-phi(i)))*Z(i);
eqDs(i)$(E0(i)<>0)..       D(i)    =e= (theta(i)**phi(i)*xid(i)*(1+tauz(i))*pz(i)
                            /pd(i))**(1/(1-phi(i)))*Z(i);
eqpzd1(i)$(E0(i)=0)..         Z(i)*(1+tauz(i))  =e=  D(i);
eqpzd2(i)$(E0(i)=0)..                    pz(i)  =e=  pd(i);
eqpzd3(i)$(E0(i)=0)..                     E(i)  =e=  0;

*[market clearing condition]
eqpqd(i)..      Q(i)    =e= sum(l,Xp(i,l)) +Xg(i) +Xv(i) +sum(j, X(i,j)) +Walras ;
*Leuis closure
eqpf(h)..       sum(j, F(h,j)) =e= FF(h);

* 3.2.6 Macro-indicators function ----------------------------------------------
eqGDP..                  GDP  =e=  sum(i,Xv(i)+Xg(i)+sum(l,Xp(i,l))+E(i)-M(i)) ;
eqnGDP..                nGDP  =e=  sum(i, pq(i)*(Xv(i)+Xg(i)+sum(l,Xp(i,l))) +pe(i)*E(i) -pm(i)*M(i) ) ;

eqGDPCHK..            GDPCHK  =e=  sum((h,j), pf(h)*F(h,j)) +sum(j, Tz(j)+Tm(j))  -sum(i, pq(i)*(Xv(i)+Xg(i)+sum(l,Xp(i,l))) +pe(i)*E(i) -pm(i)*M(i) ) ;
eqQuantityMoney..      CPI*GDP =e=  (M2+deficit*mmultiplier/2) * Velocity;
eqVelocity..          Velocity =e=  coeff_vel*sum(l,( sum(i, Xp(i,l)*pq(i)) / HOHincome(l)  ) );

eqPPI(j)..                              PPI(j)  =e=  pz(j)/ 1 *100;
eqCPI..                                    CPI  =e=  sum(j, pq(j)*sum(l,Xp0(j,l))/sum((i,l),Xp0(i,l))) / sum(j,  1 *sum(l,Xp0(j,l))/sum((i,l),Xp0(i,l)))*100;

*[fictitious objective function]
eqUU(l)..           UU(l)        =e=  sum(i, alpha(i,l)*Xp(i,l)**rhou(l))**(1/rhou(l));
eqSW..              SW           =e=  sum(l,UU(l));

eqVV(l)..           VV(l)        =e=  (HOHincome(l)  +Tp(l)  -Sp(l)  -Td(l)   -refund(l)) / sum[i,alpha(i,l)**(1/(1-rhou(l)))*pq_hoh(i,l)**(rhou(l)/(rhou(l)-1)) ]**((rhou(l)-1)/rhou(l));
eqV_source(l)..     V_source(l)  =e=  (HOHincome(l)  +Tp(l)  -Sp(l)  -Td(l)   -refund(l)) / sum[i,alpha(i,l)**(1/(1-rhou(l)))*1          **(rhou(l)/(rhou(l)-1)) ]**((rhou(l)-1)/rhou(l));
eqV_use(l)..        V_use(l)     =e=  (HOHincome0(l) +Tp0(l) -Sp0(l) -Td0(l)  -refund(l)) / sum[i,alpha(i,l)**(1/(1-rhou(l)))*pq_hoh(i,l)**(rhou(l)/(rhou(l)-1)) ]**((rhou(l)-1)/rhou(l));

eqEV(l)..           EV(l)       /100 =e=  VV(l)      *sum[i,alpha(i,l)**(1/(1-rhou(l)))*1          **(rhou(l)/(rhou(l)-1)) ]**((rhou(l)-1)/rhou(l))
                                        /(UU0(l)     *sum[i,alpha(i,l)**(1/(1-rhou(l)))*1          **(rhou(l)/(rhou(l)-1)) ]**((rhou(l)-1)/rhou(l))) -1;
eqEV_source(l)..    EV_source(l)/100 =e=  V_source(l)*sum[i,alpha(i,l)**(1/(1-rhou(l)))*1          **(rhou(l)/(rhou(l)-1)) ]**((rhou(l)-1)/rhou(l))
                                        /(UU0(l)     *sum[i,alpha(i,l)**(1/(1-rhou(l)))*1          **(rhou(l)/(rhou(l)-1)) ]**((rhou(l)-1)/rhou(l))) -1;

eqEV_use(l)..       EV_use(l)        =e=  EV(l)  -   EV_source(l);

* Assume a mild Keynesian setting: excess labor supply, but with an upward-sloping labor supply curve (finite elasticity)
eqlabordemand(hlab).. epsilonLD(hlab) * ((pf(hlab)-1)/1)    =e=   (FF(hlab)-FF0(hlab))/FF0(hlab)   ;
* ---------------------------------------------------------------------

* Initializing variables ----------------------------------------------
F.l(h,j)        =        F0(h,j)        ;
LAB.l(j)        =        LAB0(j)        ;
CAP.l(j)        =        CAP0(j)        ;
VA.l(j)        =        VA0(j)        ;
X.l(i,j)        =        X0(i,j)        ;
TX.l(j)        =        TX0(j)        ;
Z.l(j)        =        Z0(j)        ;
Xp.l(i,l)        =        Xp0(i,l)        ;
Xg.l(i)        =        Xg0(i)        ;
Xv.l(i)        =        Xv0(i)        ;
E.l(i)        =        E0(i)        ;
M.l(i)        =        M0(i)        ;
Q.l(i)        =        Q0(i)        ;
D.l(i)        =        D0(i)        ;

FF.l(h)          =       FF0(h)    ;
loan.l(l)        =        loan0(l)        ;
loang.l        =         loang0        ;
refund.l(l)        =        refund0(l)        ;
refundg.l        =          refundg0           ;
debt.l           =          debt0        ;
pF.l(h)        =        1        ;
pLAB.l(j)        =        1        ;
pCAP.l(j)        =        1        ;
pVA.l(j)        =        1        ;
pTX.l(j)        =        1        ;
pZ.l(j)        =        1        ;
pq.l(i)        =        1        ;
pq_hoh.l(i,l)      =       1       ;
pe.l(i)        =        1        ;
pm.l(i)        =        1        ;
pd.l(i)        =        1        ;
ploan.l(l)        =        1        ;
ploang.l        =        1        ;
epsilon.l        =        1        ;

HOHsalary.l(l)        =        HOHsalary0(l)        ;
HOHincome.l(l)        =        HOHincome0(l)        ;
GOVincome.l           =          GOVincome0      ;
GOVtax.l              =          GOVtax0         ;
Tp.l(l)        =        Tp0(l)        ;

Sp.l(l)        =        Sp0(l)        ;
Sg.l        =        Sg0        ;
Td.l(l)        =        Td0(l)        ;
Tz.l(j)        =        Tz0(j)        ;
Tm.l(i)        =        Tm0(i)        ;
GDP.l          =         GDP0            ;
nGDP.l          =         nGDP0            ;
GDPCHK.l       =         0               ;
PPI.l(j)        =        PPI0(j)        ;
CPI.l        =        CPI0        ;
Velocity.l        =        Velocity0        ;
Walras.l       =         0               ;

UU.l(l)        =        UU0(l)        ;
SW.l        =       SW0          ;
VV.l(l)       =       VV0(l)       ;
V_source.l(l)       =       V_source0(l)       ;
V_use.l(l)          =       V_use0(l)          ;
EV.l(l)             =       EV0(l)             ;
EV_source.l(l)      =       EV_source0(l)      ;
EV_use.l(l)         =       EV_use0(l)         ;

* ---------------------------------------------------------------------


* Numeraire and closure -----------------------------------------------
* Choose between Keynesian and Lewis-type closure. Here we fix capital endowment and allow labor supply to be elastic.
FF.fx(hcap)=FF0(hcap) ;


* ---------------------------------------------------------------------
* Defining and solving the model --------------------------------------
Model stdcge /all/;
*Solve stdcge maximizing SW using nlp;
Solve stdcge using mcp;

*$ontext
* ---------------------------------------------------------------------
set           t   scenario    /BAU,GOVexp,Investment,Corporatetaxcut,HOHtaxcut,HOHLSS,ConSubsidy/ ;
parameter
              F1(t,h,j)          "Factor demand of h by j under scenario t"
              LAB1(t,j)          "Total labor composite used by activity j under scenario t"
              CAP1(t,j)          "Total capital composite used by activity j under scenario t"
              VA1(t,j)           "Value added of activity j under scenario t"
              X1(t,i,j)          "Intermediate demand for good i by activity j under scenario t"
              TX1(t,j)           "Total intermediate input of activity j under scenario t"
              Z1(t,j)            "Gross output of activity j under scenario t"
              Xp1(t,i,l)         "Household l consumption of good i under scenario t"
              Xg1(t,i)           "Government consumption of good i under scenario t"
              Xv1(t,i)           "Investment demand for good i under scenario t"
              E1(t,i)            "Exports of good i under scenario t"
              M1(t,i)            "Imports of good i under scenario t"
              Q1(t,i)            "Armington composite quantity of good i under scenario t"
              D1(t,i)            "Domestic supply of good i under scenario t"

              pF1(t,h)           "Factor price of h under scenario t"
              pLAB1(t,j)         "Price of labor composite in activity j under scenario t"
              pCAP1(t,j)         "Price of capital composite in activity j under scenario t"
              pVA1(t,j)          "Price of value added in activity j under scenario t"
              pTX1(t,j)          "Price index of intermediates for activity j under scenario t"
              pZ1(t,j)           "Producer price of activity j under scenario t"
              pq1(t,i)           "Price of Armington composite good i under scenario t"
              pq_hoh1(t,i,l)     "Consumer price of good i for household l under scenario t"
              pe1(t,i)           "Export price of good i in domestic currency under scenario t"
              pm1(t,i)           "Import price of good i in domestic currency under scenario t"
              pd1(t,i)           "Domestic price of good i under scenario t"
              epsilon1(t)        "Exchange rate under scenario t"

              FF1(t,h)           "Factor endowment/supply of factor h under scenario t"
              ploan1(t,l)        "Loan price (1+interest) to household l under scenario t"
              loan1(t,l)         "New borrowing of household l under scenario t"
              refund1(t,l)       "Debt service payment of household l under scenario t"
              loang1(t)          "Government borrowing under scenario t"
              refundg1(t)        "Government debt service payment under scenario t"

              GOVincome1(t)      "Total government income under scenario t"
              GOVtax1(t)         "Total tax revenue under scenario t"

              HOHincome1(t,l)    "Total income of household l under scenario t"
              HOHsalary1(t,l)    "Factor income of household l under scenario t"
              Tp1(t,l)           "Government transfers to household l under scenario t"

              Sp1(t,l)           "Private saving of household l under scenario t"
              Sg1(t)             "Government saving under scenario t"
              Td1(t,l)           "Direct tax paid by household l under scenario t"
              Tz1(t,j)           "Production tax revenue from activity j under scenario t"
              Tm1(t,i)           "Import tariff revenue on good i under scenario t"

              GDP1(t)            "Real GDP under scenario t"
              nGDP1(t)           "Nominal GDP under scenario t"
              GDPCHK1(t)         "GDP accounting check under scenario t"
              PPI1(t,j)          "Producer price index for activity j under scenario t"
              CPI1(t)            "Consumer price index under scenario t"
              Velocity1(t)       "Velocity of money under scenario t"
              Walras1(t)         "Walrasian slack under scenario t"

              UU1(t,l)           "Utility of household l under scenario t"
              SW1(t)             "Social welfare under scenario t"
              VV1(t,l)           "Money-metric utility of household l under scenario t"
              V_source1(t,l)     "Source-side utility index of household l under scenario t"
              V_use1(t,l)        "Use-side utility index of household l under scenario t"
              EV1(t,l)           "Equivalent variation of household l under scenario t"
              EV_source1(t,l)    "Source-side EV of household l under scenario t"
              EV_use1(t,l)       "Use-side EV of household l under scenario t"
;

Loop(t,
F1(t,h,j)        =        F.l(h,j)        ;
LAB1(t,j)        =        LAB.l(j)        ;
CAP1(t,j)        =        CAP.l(j)        ;
VA1(t,j)        =        VA.l(j)        ;
X1(t,i,j)        =        X.l(i,j)        ;
TX1(t,j)        =        TX.l(j)        ;
Z1(t,j)        =        Z.l(j)        ;
Xp1(t,i,l)        =        Xp.l(i,l)        ;
Xg1(t,i)        =        Xg.l(i)        ;
Xv1(t,i)        =        Xv.l(i)        ;
E1(t,i)        =        E.l(i)        ;
M1(t,i)        =        M.l(i)        ;
Q1(t,i)        =        Q.l(i)        ;
D1(t,i)        =        D.l(i)        ;

FF1(t,h)       =         FF.l(h)       ;
pF1(t,h)        =        pF.l(h)        ;
pLAB1(t,j)        =        pLAB.l(j)        ;
pCAP1(t,j)        =        pCAP.l(j)        ;
pVA1(t,j)        =        pVA.l(j)        ;
pTX1(t,j)        =        pTX.l(j)        ;
pZ1(t,j)        =        pZ.l(j)        ;
pq1(t,i)        =        pq.l(i)        ;
pq_hoh1(t,i,l)        =        pq_hoh.l(i,l)        ;
pe1(t,i)        =        pe.l(i)        ;
pm1(t,i)        =        pm.l(i)        ;
pd1(t,i)        =        pd.l(i)        ;
epsilon1(t)        =        epsilon.l        ;

ploan1(t,l)        =        ploan.l(l)        ;
loan1(t,l)        =        loan.l(l)        ;
refund1(t,l)      =        refund.l(l)      ;
loang1(t)         =        loang.l       ;
refundg1(t)     =        refundg.l     ;

GOVincome1(t)     =        GOVincome.l     ;
GOVtax1(t)        =        GOVtax.l        ;

HOHincome1(t,l)        =        HOHincome.l(l)        ;
HOHsalary1(t,l)        =         HOHsalary.l(l)       ;
Tp1(t,l)        =        Tp.l(l)        ;

Sp1(t,l)        =        Sp.l(l)        ;
Sg1(t)        =        Sg.l        ;
Td1(t,l)        =        Td.l(l)        ;
Tz1(t,j)        =        Tz.l(j)        ;
Tm1(t,i)        =        Tm.l(i)        ;

GDP1(t)        =        GDP.l        ;
nGDP1(t)        =        nGDP.l        ;
GDPCHK1(t)        =        GDPCHK.l        ;
PPI1(t,j)        =        PPI.l(j)        ;
CPI1(t)        =        CPI.l        ;
Velocity1(t)        =        Velocity.l        ;
Walras1(t)        =        Walras.l        ;

UU1(t,l)        =        UU.l(l)        ;
SW1(t)        =        SW.l        ;
VV1(t,l)       =     VV.l(l)       ;
V_source1(t,l)       =     V_source.l(l)       ;
V_use1(t,l)          =     V_use.l(l)          ;
EV1(t,l)             =     EV.l(l)             ;
EV_source1(t,l)      =     EV_source.l(l)      ;
EV_use1(t,l)         =     EV_use.l(l)         ;

* Policy shock: government issues 1 trillion RMB of new bonds (deficit = 1000)--
deficit =       1000;

if(      ord(t)  =       1,
         Solve stdcge using mcp;
       );

if(      ord(t)  =       2,
         investment=     deficit;
         Solve stdcge using mcp;
         investment=     0;
       );

if(      ord(t)  =       3,
         tauz(i)  =       tauz(i)*(1-deficit/sum(j,Tz0(j)));
         Solve stdcge using mcp;
         tauz(i)  =       tauz(i)/(1-deficit/sum(j,Tz0(j)));
       );

if(      ord(t)  =       4,
         taud(l) =       taud(l)*(1-deficit/sum(al,Td0(al)));
         Solve stdcge using mcp;
         taud(l) =       taud(l)/(1-deficit/sum(al,Td0(al)));
       );

if(      ord(t)  =       5,
         LSS(l)   =       deficit/card(l);
         Solve stdcge using mcp;
         LSS(l)   =       0;
       );

if(      ord(t)  =       6,
         coupon  =       deficit;
         share_coupon(i,l)= 0;
         share_coupon("PFEDN",l)= 0.8/card(l) ;
         share_coupon("HOU",l)= 0.2/card(l) ;
         Solve stdcge using mcp;
         coupon  =       0;
         share_coupon(i,l)= Xp0(i,l)/sum((j,al), Xp0(j,al));
       );
);

Parameter
               dFF1(t,h)           "Percentage change in factor endowment of h vs BAU under t"
               dF1(t,h,j)          "Percentage change in factor demand F(h,j) vs BAU under t"
               dLAB1(t,j)          "Percentage change in labor composite LAB(j) vs BAU under t"
               dCAP1(t,j)          "Percentage change in capital composite CAP(j) vs BAU under t"
               dVA1(t,j)           "Percentage change in value added VA(j) vs BAU under t"
               dX1(t,i,j)          "Percentage change in intermediate input X(i,j) vs BAU under t"
               dTX1(t,j)           "Percentage change in TX(j) vs BAU under t"
               dZ1(t,j)            "Percentage change in output Z(j) vs BAU under t"
               dXp1(t,i,l)         "Percentage change in household consumption Xp(i,l) vs BAU under t"
               dXg1(t,i)           "Percentage change in government demand Xg(i) vs BAU under t"
               dXv1(t,i)           "Percentage change in investment demand Xv(i) vs BAU under t"
               dE1(t,i)            "Percentage change in exports E(i) vs BAU under t"
               dM1(t,i)            "Percentage change in imports M(i) vs BAU under t"
               dQ1(t,i)            "Percentage change in Armington composite Q(i) vs BAU under t"
               dD1(t,i)            "Percentage change in domestic supply D(i) vs BAU under t"

               dFF1(t,h)           "Percentage change in factor supply FF(h) vs BAU under t"
               dpF1(t,h)           "Percentage change in factor price pF(h) vs BAU under t"
               dpLAB1(t,j)         "Percentage change in pLAB(j) vs BAU under t"
               dpCAP1(t,j)         "Percentage change in pCAP(j) vs BAU under t"
               dpVA1(t,j)          "Percentage change in pVA(j) vs BAU under t"
               dpTX1(t,j)          "Percentage change in pTX(j) vs BAU under t"
               dpZ1(t,j)           "Percentage change in pZ(j) vs BAU under t"
               dpq1(t,i)           "Percentage change in pq(i) vs BAU under t"
               dpq_hoh1(t,i,l)     "Percentage change in pq_hoh(i,l) vs BAU under t"
               dpe1(t,i)           "Percentage change in export price pe(i) vs BAU under t"
               dpm1(t,i)           "Percentage change in import price pm(i) vs BAU under t"
               dpd1(t,i)           "Percentage change in domestic price pd(i) vs BAU under t"
               depsilon1(t)        "Percentage change in exchange rate epsilon vs BAU under t"

               dploan1(t,l)        "Percentage change in loan price to household l vs BAU"
               dloan1(t,l)         "Percentage change in borrowing of household l vs BAU"
               drefund1(t,l)       "Percentage change in debt service of household l vs BAU"
               dloang1(t)          "Percentage change in government borrowing vs BAU"
               drefundg1(t)        "Percentage change in government debt service vs BAU"

               dGOVincome1(t)      "Percentage change in government income vs BAU under t"
               dGOVtax1(t)         "Percentage change in government tax revenue vs BAU under t"

               dHOHincome1(t,l)    "Percentage change in HOHincome(l) vs BAU under t"
               dHOHsalary1(t,l)    "Percentage change in HOHsalary(l) vs BAU under t"
               dTp1(t,l)           "Percentage change in transfers Tp(l) vs BAU under t"

               dSp1(t,l)           "Percentage change in private saving Sp(l) vs BAU under t"
               dSg1(t)             "Percentage change in government saving vs BAU under t"
               dTd1(t,l)           "Percentage change in Td(l) vs BAU under t"
               dTz1(t,j)           "Percentage change in Tz(j) vs BAU under t"
               dTm1(t,i)           "Percentage change in Tm(i) vs BAU under t"

               dGDP1(t)            "Percentage change in real GDP vs BAU under t"
               dnGDP1(t)           "Percentage change in nominal GDP vs BAU under t"
               dGDPCHK1(t)         "Change in GDPCHK vs BAU (accounting check)"
               dPPI1(t,i)          "Percentage change in PPI(i) vs BAU under t"
               dCPI1(t)            "Percentage change in CPI vs BAU under t"
               dVelocity1(t)       "Percentage change in money velocity vs BAU under t"
               dWalras1(t)         "Change in Walrasian slack vs BAU under t"

               dUU1(t,l)           "Percentage change in utility UU(l) vs BAU under t"
               dSW1(t)             "Percentage change in social welfare SW vs BAU under t"
               dVV1(t,l)           "Percentage change in VV(l) vs BAU under t"
               dV_source1(t,l)     "Percentage change in V_source(l) vs BAU under t"
               dV_use1(t,l)        "Percentage change in V_use(l) vs BAU under t"
               dEV1(t,l)           "Change in EV(l) vs BAU (if not normalized)"
               dEV_source1(t,l)    "Change in EV_source(l) vs BAU"
               dEV_use1(t,l)       "Change in EV_use(l) vs BAU"
;

dFF1(t,h)         =(FF1(t,h)/FF1("BAU",h)-1)*100;
dF1(t,h,j)        =(F1(t,h,j)/F1("BAU",h,j)-1)*100;
dLAB1(t,j)        =(LAB1(t,j)/LAB1("BAU",j)-1)*100;
dCAP1(t,j)        =(CAP1(t,j)/CAP1("BAU",j)-1)*100;
dVA1(t,j)        =(VA1(t,j)/VA1("BAU",j)-1)*100;
dX1(t,i,j)        =(X1(t,i,j)/X1("BAU",i,j)-1)*100;
dTX1(t,j)        =(TX1(t,j)/TX1("BAU",j)-1)*100;
dZ1(t,j)        =(Z1(t,j)/Z1("BAU",j)-1)*100;
dXp1(t,i,l)        =(Xp1(t,i,l)/Xp1("BAU",i,l)-1)*100;
dXg1(t,i)$Xg0(i) =(Xg1(t,i)/Xg1("BAU",i)-1)*100;
dXv1(t,i)$Xv0(i) =(Xv1(t,i)/Xv1("BAU",i)-1)*100;
dE1(t,i)$E0(i)  =(E1(t,i)/E1("BAU",i)-1)*100;
dM1(t,i)$M0(i)  =(M1(t,i)/M1("BAU",i)-1)*100;
dQ1(t,i)        =(Q1(t,i)/Q1("BAU",i)-1)*100;
dD1(t,i)        =(D1(t,i)/D1("BAU",i)-1)*100;

dFF1(t,h)        =(FF1(t,h)/FF1("BAU",h)-1)*100;
dpF1(t,h)        =(pF1(t,h)/pF1("BAU",h)-1)*100;
dpLAB1(t,j)        =(pLAB1(t,j)/pLAB1("BAU",j)-1)*100;
dpCAP1(t,j)        =(pCAP1(t,j)/pCAP1("BAU",j)-1)*100;
dpVA1(t,j)        =(pVA1(t,j)/pVA1("BAU",j)-1)*100;
dpTX1(t,j)        =(pTX1(t,j)/pTX1("BAU",j)-1)*100;
dpZ1(t,j)        =(pZ1(t,j)/pZ1("BAU",j)-1)*100;
dpq1(t,i)        =(pq1(t,i)/pq1("BAU",i)-1)*100;
dpq_hoh1(t,i,l)        =(pq_hoh1(t,i,l)/pq_hoh1("BAU",i,l)-1)*100;
dpe1(t,i)        =(pe1(t,i)/pe1("BAU",i)-1)*100;
dpm1(t,i)        =(pm1(t,i)/pm1("BAU",i)-1)*100;
dpd1(t,i)        =(pd1(t,i)/pd1("BAU",i)-1)*100;
depsilon1(t)        =(epsilon1(t)/epsilon1("BAU")-1)*100;

dploan1(t,l)      =   (ploan1(t,l)    / ploan1("BAU",l)   -1)*100;
dloan1(t,l)       =   (loan1(t,l)     / loan1("BAU",l)    -1)*100;
drefund1(t,l)     =   (refund1(t,l)   / refund1("BAU",l)  -1)*100;
dloang1(t)        =   (loang1(t)      / loang1("BAU")     -1)*100;
drefundg1(t)      =   (refundg1(t)    / refundg1("BAU")   -1)*100;

dGOVincome1(t)    =   (GOVincome1(t)  / GOVincome1("BAU") -1)*100;
dGOVtax1(t)       =   (GOVtax1(t)     / GOVtax1("BAU")    -1)*100;

dHOHincome1(t,l)        =(HOHincome1(t,l)/HOHincome1("BAU",l)-1)*100;
dHOHsalary1(t,l)        =(HOHsalary1(t,l)/HOHsalary1("BAU",l)-1)*100;
dTp1(t,l)        =(Tp1(t,l)/Tp1("BAU",l)-1)*100;

dSp1(t,l)        =(Sp1(t,l)/Sp1("BAU",l)-1)*100;
dSg1(t)        =(Sg1(t)/Sg1("BAU")-1)*100;
dTd1(t,l)$Td0(l) =(Td1(t,l)/Td1("BAU",l)-1)*100;
dTz1(t,j)        =(Tz1(t,j)/Tz1("BAU",j)-1)*100;
dTm1(t,i)$Tm0(i) =(Tm1(t,i)/Tm1("BAU",i)-1)*100;

dGDP1(t)        =(GDP1(t)/GDP1("BAU")-1)*100;
dnGDP1(t)       =(nGDP1(t)/nGDP1("BAU")-1)*100;
dPPI1(t,j)        =(PPI1(t,j)/PPI1("BAU",j)-1)*100;
dCPI1(t)        =(CPI1(t)/CPI1("BAU")-1)*100;
dVelocity1(t)        =(Velocity1(t)/Velocity1("BAU")-1)*100;

dUU1(t,l)                        =(UU1(t,l)/UU1("BAU",l)-1)*100;
dSW1(t)                          =(SW1(t)/SW1("BAU")-1)*100;

execute_unload 'output.gdx';
