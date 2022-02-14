options
profile = 1
limrow = 10000
limcol = 10000
optCR = 1E-9
optCA = 0.0
iterlim = 100000000
decimals = 8
reslim = 1000000
work = 500000000
sysout = off
solprint = on
lp = cplex
mip = cplex
;

******************************** SET DEFINITIONS ***************************************
*
*       i                       = set of all metabolites
*       j                       = set of all reactions
*
*****************************************************************************************

$setnames "%gams.i%" filepath filename fileextension

Sets

i
/
$include %filepath%gams_inputs/metabolites.txt
/

j
/
$include %filepath%gams_inputs/reactions.txt

/

medium(j)
/
$include %filepath%gams_inputs/medium.txt
/

exchange(j)
/
$include %filepath%gams_inputs/ExchangeRxns.txt
/

rxns_metabcons
/
$include %filepath%gams_inputs/RxnsWithMetabCons_AllMets.txt
/

rxns_SNPcons
/
$include %filepath%gams_inputs/RxnsWithSNPs.txt
/

l
/
$include %filepath%gams_inputs/AllGenotypesInStudy.txt
/

k
/
$include %filepath%gams_inputs/SNPList_AllEcotypes.txt
/

dummy(j)

;


****************************** PARAMETRIC CONSTANTS USED ********************************
*
*       s(i,j)          = Stoichiometry of metabolite i in reaction j
*       LB(j)/UB(j)     = Stores the lower and upper bounds for each reaction j
*       LB_pfba(j)/UB_pfba(j)     = Stores the lower and upper bounds for each reaction j from pFBA at max biomass
*		Snp2EcoMap(k,l) = maps SNP k in genotype l (1 - present, 0 - absent)
*		Rxn2Snp2EcoMap = maps if reaction j in genotype l has a SNP k in associated genes (1 - present, 0 - absent)
*
*****************************************************************************************
Parameters

S(i,j)
/
$include %filepath%gams_inputs/sij_matrix.txt
/

LB(j) , UB(j)
LB_pfba(j) , UB_pfba(j)

Snp2EcoMap(k,l)
Rxn2Snp2EcoMap(j,k,l)

model_stat

;

***************************** VARIABLE DEFINITIONS **************************************
*
*		v_ref(j)        = Flux of a reaction j (+/-) in the reference genotype
*       v(j,l)          = Flux of a reaction j (+/-) in genotype l
*       z               = Objective value (+/-)
*		S_pos(j,l)		= slack due to SNPs (LB_pfba limiting)
*		S_neg(j,l)		= slack due to SNPs (UB_pfba limiting)
*		MassActSlack	= slack due to deviations from imposed mass action kinetics (+/-)
*		MassActSlack_pos(j,l)	= slack due to deviations from imposed mass action kinetics (+)
*		max_bio			= maximum biomass for the reference genotype
*		X_pos(k,l)      = flux change due to SNP k in genotype l (inactivating SNP)
*		X_neg(k,l)      = flux change due to SNP k in genotype l (activating SNP)
*		y_pos(k) 		= binary that is 1 when SNP k is inactivating, 0 otherwise
*		y_neg(k) 		= binary that is 1 when SNP k is inactivating, 0 otherwise
*
*****************************************************************************************


Variables

z
v(j,l)
v_ref(j)

S_pos(j,l)
S_neg(j,l)

MassActSlack(j,l)
MassActSlack_pos(j,l)

tot_min_MassActslack

max_bio
tot_min_slack
X_pos(k,l)
X_neg(k,l)

;

Positive Variables

S_pos(j,l)
S_neg(j,l)
X_pos(k,l)
X_neg(k,l)
;


Binary Variables

y_pos(k)
y_neg(k)

;

***************************** EQUATION DEFINITIONS **************************************
*
*       Stoic(i)        = Stoichiometric Constraint
*       Obj             = Maximize/Minimize flux of each reaction one at a time

*****************************************************************************************
Equations

Obj_maxBio
Obj_minMassActSlack
Stoic(i,l)
Stoic_ref(i)

$include %filepath%gams_inputs/Biomass_ConsList_AllGenotypes.txt
$include %filepath%gams_inputs/MetabConsList_MinMaxMets_AllEcotypes_BwdRxnsflipped.txt
$include %filepath%gams_inputs/MetabConsList_measuredMets_AllEcotypes.txt
$include %filepath%gams_inputs/RxnBndsList_MetSNPseparated.txt

carbonuptake_auto
carbonuptake_auto_allecos(l)

DecomposePosSlack(j,l)
DecomposeNegSlack(j,l)

ConsSNPBehaviour_pos_ub(k,l)
ConsSNPBehaviour_pos_lb(k,l)
ConsSNPBehaviour_neg_ub(k,l)
ConsSNPBehaviour_neg_lb(k,l)
BinSumCons(k)

MassActSlack_pos_ub(j,l)
MassActSlack_pos_lb(j,l)

tot_min_MassActslackCons

;

****add the mapping matrices 
Snp2EcoMap(k,l) = 0;
$include %filepath%gams_inputs/Snp2EcoMapping.txt

Rxn2Snp2EcoMap(j,k,l) = 0;
$include %filepath%gams_inputs/Rxn2Snp2EcoMapping.txt

****DEFININING THE CONSTRAINTS *************************************************************

Obj_minMassActSlack..           z =e= sum(l, sum(j, 100*MassActSlack_pos(j,l) +  S_pos(j,l) + S_neg(j,l) ));
Obj_maxBio..					z =e= v_ref('BiomassRxn_Leaf');

Stoic(i,l)..					sum(j, S(i,j) * v(j,l) )  =e= 0;
Stoic_ref(i)..					sum(j, S(i,j) * v_ref(j) )  =e= 0;

carbonuptake_auto..				v_ref('EX_carbon_dioxide')  =g= -10.0;
carbonuptake_auto_allecos(l)..	v('BiomassRxn_Leaf',l) =g= -10;

$include %filepath%gams_inputs/Biomass_Cons_AllGenotypes.txt

***with mass action slack 
$include %filepath%gams_inputs/MetabCons_MinMaxMets_AllEcotypes_BwdExnsflipped.txt
$include %filepath%gams_inputs/MetabCons_measuredMets_AllEcotypes.txt
$include %filepath%gams_inputs/RxnBnds_MetSNPseparated.txt

DecomposePosSlack(j,l)..							S_pos(j,l) =e= sum(k, Rxn2Snp2EcoMap(j,k,l)*X_pos(k,l));
DecomposeNegSlack(j,l)..							S_neg(j,l) =e= sum(k, Rxn2Snp2EcoMap(j,k,l)*X_neg(k,l));

ConsSNPBehaviour_pos_ub(k,l)..						X_pos(k,l) =l= (100.0)*Snp2EcoMap(k,l)*y_pos(k);
ConsSNPBehaviour_pos_lb(k,l)..						X_pos(k,l) =g= (0.0)*Snp2EcoMap(k,l)*y_pos(k) ;
ConsSNPBehaviour_neg_ub(k,l)..						X_neg(k,l) =l= (100.0)*Snp2EcoMap(k,l)*y_neg(k) ;
ConsSNPBehaviour_neg_lb(k,l)..						X_neg(k,l) =g= (0.0)*Snp2EcoMap(k,l)*y_neg(k) ;

BinSumCons(k)..										y_pos(k) + y_neg(k) =l= 1;

MassActSlack_pos_ub(j,l)..							MassActSlack_pos(j,l) =g= MassActSlack(j,l);
MassActSlack_pos_lb(j,l)..							MassActSlack_pos(j,l) =g= -MassActSlack(j,l);

tot_min_MassActslackCons..		sum(l, sum(j, MassActSlack_pos(j,l))) =l= tot_min_MassActslack;

*****************************************************************************************
***** DEFINING THE UPPER AND LOWER BOUNDS FOR ALL THE REACTIONS *************************


****the reference flux distribution 
***populates LB_pfba(j) , UB_pfba(j)
$include %filepath%gams_inputs/reaction_bounds_pFBA.txt

***load bounds for all other vars  
$include %filepath%gams_inputs/reaction_bounds.txt

UB(j)$(exchange(j)) = 10000;
LB(j)$(exchange(j)) = 0;

LB(j)$(medium(j)) = -10000;

LB('ATPM_c') = 0.6772;
UB('ATPM_c') = 0.6772;

LB('BiomassRxn_Leaf') = 0;
UB('BiomassRxn_Leaf') = 10000;

LB('EX_carbon_dioxide') = -10000;

v_ref.lo(j)$exchange(j) = 0.0;
v_ref.up(j)$exchange(j) = 10000.0;
***medium
v_ref.lo(j)$medium(j) = -1000.0;
***add co2
v_ref.lo('EX_carbon_dioxide') = -10;

*** SETTING LOWER AND UPPER BOUNDS FOR THE FLUXES (variable bounds)
v.lo(j,l) = LB(j);
v.up(j,l) = UB(j);

v_ref.lo(j) = LB(j);
v_ref.up(j) = UB(j);

v.lo(j,l)$exchange(j) = 0.0;
v.up(j,l)$exchange(j) = 10000.0;
***medium
v.lo(j,l)$medium(j) = -10000.0;



************************ DECLARING THE MODEL with constraints****************************

Model fba_auto
/
Obj_maxBio
Stoic_ref
carbonuptake_auto
/;
fba_auto.optfile = 1;

Model minMassActSlack
/
Obj_minMassActSlack
Stoic
$include %filepath%gams_inputs/MetabConsList_measuredMets_AllEcotypes.txt
$include %filepath%gams_inputs/MetabConsList_MinMaxMets_AllEcotypes_BwdRxnsflipped.txt
$include %filepath%gams_inputs/Biomass_ConsList_AllGenotypes.txt
$include %filepath%gams_inputs/RxnBndsList_MetSNPseparated.txt
DecomposePosSlack
DecomposeNegSlack
ConsSNPBehaviour_pos_ub
ConsSNPBehaviour_pos_lb
ConsSNPBehaviour_neg_ub
ConsSNPBehaviour_neg_lb
BinSumCons
MassActSlack_pos_ub
MassActSlack_pos_lb
/;
minMassActSlack.optfile = 1;
minMassActSlack.holdfixed = 1;



file report /results.txt/;

put report

alias(j,j1);
alias(j,j2);

***make light not limiting
v.lo('EX_light',l) = -10000;
v.up('EX_proton',l) = 10000;
v_ref.lo('EX_light') = -10000;
v_ref.up('EX_proton') = 10000;

***slack variable bounds
S_pos.lo(j,l) = 0;
S_pos.up(j,l) = 100;
S_neg.lo(j,l) = 0;
S_neg.up(j,l) = 100;

X_pos.lo(k,l) = 0;
X_pos.up(k,l) = 100;
X_neg.lo(k,l) = 0;
X_neg.up(k,l) = 100;

MassActSlack.lo(j,l) = -1000;
MassActSlack.up(j,l) = 1000;

***************Solve the FBA first, to get max biomass in the reference genotype***************

Solve fba_auto using lp maximizing z;
model_stat = fba_auto.modelstat;
max_bio.fx = z.l;
put "objective flux:","        ",z.l:10:5,"    ",model_stat/;       
put "biomass prod (control or auto):",v_ref.l('BiomassRxn_Leaf'):10:5/;
put "flux_r2(co2):   ",v_ref.l('EX_carbon_dioxide'):10:5/;



S_pos.fx(j,l)$(not rxns_SNPcons(j)) = 0;
S_neg.fx(j,l)$(not rxns_SNPcons(j)) = 0;
S_pos.fx('BiomassRxn_Leaf',l) = 0.0;
S_neg.fx('BiomassRxn_Leaf',l) = 0.0;

Solve minMassActSlack using mip minimizing z;
model_stat = minMassActSlack.modelstat;
put "objective flux (tot min mass action slack):","        ",z.l:10:5,"    ",model_stat/;
put "tot min SNPdev:","        ",sum(l,sum(j, (S_pos.l(j, l) + S_neg.l(j, l))))/;
put "tot min MassActSlack:","        ",sum(l, sum(j, MassActSlack_pos.l(j,l)))/;
tot_min_MassActslack.fx  = z.l;


****print results now!
loop(l,
put l.tl:40,", Biomass: ",v.l('BiomassRxn_Leaf',l):10:5/;
);

loop(l,
put //"Looking in Ecotype: ",l.tl:0/;

loop(j$(S_pos.l(j,l) ne 0),
put j.tl:0,", S_pos: ",S_pos.l(j,l):0:10/;);

loop(j$(S_neg.l(j,l) ne 0),
put j.tl:0,", S_neg: ",S_neg.l(j,l):0:10/;);

put //;

loop(k$(X_pos.l(k,l) ne 0.0),
put k.tl:0,"     , X_pos: ",X_pos.l(k,l):0:10,"    ",y_pos.l(k):10:5/;);

loop(k$(X_neg.l(k,l) ne 0),
put k.tl:0,"     , X_neg: ",X_neg.l(k,l):0:10,"    ",y_neg.l(k):10:5/;
);

put //"Exchange rxns now----"//;
loop(j$(exchange(j)),
put j.tl:0,"   ",v.l(j,l):0:10/;);

);

put //;

