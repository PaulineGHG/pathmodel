%*

This is work in progress on the oxylipin pathway in brown algae. Initial code by Enora Corre.

Data
Describing the starting set of metabolites

atom(Compound name, Atom number, Atom type).
bond(Compound name, Chemical bond type, Atom number1, Greater Atom number2).

Warning: in the atom description, "glut" stands provisionally for glutathion. Probably will need to be properly encoded as a domain.

*%

% linoleic_acid
atom("linoleic_acid",1..18,carb). atom("linoleic_acid",19..20,oxyg).
bond("linoleic_acid",single,1,2). bond("linoleic_acid",single,2,3). bond("linoleic_acid",single,3,4).
bond("linoleic_acid",single,4,5). bond("linoleic_acid",single,5,6). bond("linoleic_acid",single,6,7).
bond("linoleic_acid",single,7,8). bond("linoleic_acid",single,8,9). bond("linoleic_acid",double,9,10).
bond("linoleic_acid",single,10,11). bond("linoleic_acid",single,11,12). bond("linoleic_acid",double,12,13).
bond("linoleic_acid",single,13,14). bond("linoleic_acid",single,14,15). bond("linoleic_acid",single,15,16).
bond("linoleic_acid",single,16,17). bond("linoleic_acid",single,17,18). bond("linoleic_acid",single,1,19).
bond("linoleic_acid",double,1,20).

% 13-HpODE
atom("13-HpODE",1..18,carb). atom("13-HpODE",19..22,oxyg).
bond("13-HpODE",single,1,2). bond("13-HpODE",single,2,3). bond("13-HpODE",single,3,4).
bond("13-HpODE",single,4,5). bond("13-HpODE",single,5,6). bond("13-HpODE",single,6,7).
bond("13-HpODE",single,7,8). bond("13-HpODE",single,8,9). bond("13-HpODE",double,9,10).
bond("13-HpODE",single,10,11). bond("13-HpODE",double,11,12). bond("13-HpODE",single,12,13).
bond("13-HpODE",single,13,14). bond("13-HpODE",single,14,15). bond("13-HpODE",single,15,16).
bond("13-HpODE",single,16,17). bond("13-HpODE",single,17,18). bond("13-HpODE",single,1,19).
bond("13-HpODE",double,1,20). bond("13-HpODE",single,13,21). bond("13-HpODE",single,21,22).

% 12,13-EpOME
atom("12,13-EpOME",1..18,carb). atom("12,13-EpOME",19..21,oxyg).
bond("12,13-EpOME",single,1,2). bond("12,13-EpOME",single,2,3). bond("12,13-EpOME",single,3,4).
bond("12,13-EpOME",single,4,5). bond("12,13-EpOME",single,5,6). bond("12,13-EpOME",single,6,7).
bond("12,13-EpOME",single,7,8). bond("12,13-EpOME",single,8,9). bond("12,13-EpOME",double,9,10).
bond("12,13-EpOME",single,10,11). bond("12,13-EpOME",single,11,12). bond("12,13-EpOME",single,12,13).
bond("12,13-EpOME",single,13,14). bond("12,13-EpOME",single,14,15). bond("12,13-EpOME",single,15,16).
bond("12,13-EpOME",single,16,17). bond("12,13-EpOME",single,17,18). bond("12,13-EpOME",single,1,19).
bond("12,13-EpOME",double,1,20). bond("12,13-EpOME",single,12,21). bond("12,13-EpOME",single,13,21).

% 12-hydroxy-13-glutathion-OME
atom("12-hydroxy-13-glutathion-OME",1..18,carb). atom("12-hydroxy-13-glutathion-OME",19..21,oxyg). atom("12-hydroxy-13-glutathion-OME",22,glut).
bond("12-hydroxy-13-glutathion-OME",single,1,2). bond("12-hydroxy-13-glutathion-OME",single,2,3). bond("12-hydroxy-13-glutathion-OME",single,3,4).
bond("12-hydroxy-13-glutathion-OME",single,4,5). bond("12-hydroxy-13-glutathion-OME",single,5,6). bond("12-hydroxy-13-glutathion-OME",single,6,7).
bond("12-hydroxy-13-glutathion-OME",single,7,8). bond("12-hydroxy-13-glutathion-OME",single,8,9). bond("12-hydroxy-13-glutathion-OME",double,9,10).
bond("12-hydroxy-13-glutathion-OME",single,10,11). bond("12-hydroxy-13-glutathion-OME",single,11,12). bond("12-hydroxy-13-glutathion-OME",single,12,13).
bond("12-hydroxy-13-glutathion-OME",single,13,14). bond("12-hydroxy-13-glutathion-OME",single,14,15). bond("12-hydroxy-13-glutathion-OME",single,15,16).
bond("12-hydroxy-13-glutathion-OME",single,16,17). bond("12-hydroxy-13-glutathion-OME",single,17,18). bond("12-hydroxy-13-glutathion-OME",single,1,19).
bond("12-hydroxy-13-glutathion-OME",double,1,20). bond("12-hydroxy-13-glutathion-OME",single,12,21). bond("12-hydroxy-13-glutathion-OME",single,13,22).

% 9-HPOD
atom("9-HPOD",1..18,carb). atom("9-HPOD",19..22,oxyg).
bond("9-HPOD",single,1,2). bond("9-HPOD",single,2,3). bond("9-HPOD",single,3,4).
bond("9-HPOD",single,4,5). bond("9-HPOD",single,5,6). bond("9-HPOD",single,6,7).
bond("9-HPOD",single,7,8). bond("9-HPOD",single,8,9). bond("9-HPOD",single,9,10).
bond("9-HPOD",double,10,11). bond("9-HPOD",single,11,12). bond("9-HPOD",double,12,13).
bond("9-HPOD",single,13,14). bond("9-HPOD",single,14,15). bond("9-HPOD",single,15,16).
bond("9-HPOD",single,16,17). bond("9-HPOD",single,17,18). bond("9-HPOD",single,1,19).
bond("9-HPOD",double,1,20). bond("9-HPOD",single,9,21). bond("9-HPOD",single,21,22).

% 9,10-epoxy-11-hydroxy-12-octdecenoic 
atom("9,10-epoxy-11-hydroxy-12-octdecenoic",1..18,carb). atom("9,10-epoxy-11-hydroxy-12-octdecenoic",19..22,oxyg).
bond("9,10-epoxy-11-hydroxy-12-octdecenoic",single,1,2). bond("9,10-epoxy-11-hydroxy-12-octdecenoic",single,2,3). bond("9,10-epoxy-11-hydroxy-12-octdecenoic",single,3,4).
bond("9,10-epoxy-11-hydroxy-12-octdecenoic",single,4,5). bond("9,10-epoxy-11-hydroxy-12-octdecenoic",single,5,6). bond("9,10-epoxy-11-hydroxy-12-octdecenoic",single,6,7).
bond("9,10-epoxy-11-hydroxy-12-octdecenoic",single,7,8). bond("9,10-epoxy-11-hydroxy-12-octdecenoic",single,8,9). bond("9,10-epoxy-11-hydroxy-12-octdecenoic",single,9,10).
bond("9,10-epoxy-11-hydroxy-12-octdecenoic",single,10,11). bond("9,10-epoxy-11-hydroxy-12-octdecenoic",single,11,12). bond("9,10-epoxy-11-hydroxy-12-octdecenoic",double,12,13).
bond("9,10-epoxy-11-hydroxy-12-octdecenoic",single,13,14). bond("9,10-epoxy-11-hydroxy-12-octdecenoic",single,14,15). bond("9,10-epoxy-11-hydroxy-12-octdecenoic",single,15,16).
bond("9,10-epoxy-11-hydroxy-12-octdecenoic",single,16,17). bond("9,10-epoxy-11-hydroxy-12-octdecenoic",single,17,18). bond("9,10-epoxy-11-hydroxy-12-octdecenoic",single,1,19).
bond("9,10-epoxy-11-hydroxy-12-octdecenoic",double,1,20). bond("9,10-epoxy-11-hydroxy-12-octdecenoic",single,9,21). bond("9,10-epoxy-11-hydroxy-12-octdecenoic",single,10,21).
bond("9,10-epoxy-11-hydroxy-12-octdecenoic",single,11,22).

% leukotrieneA4
atom("leukotrieneA4",1..20,carb). atom("leukotrieneA4",21..23,oxyg).
bond("leukotrieneA4",single,1,2). bond("leukotrieneA4",single,2,3). bond("leukotrieneA4",single,3,4).
bond("leukotrieneA4",single,4,5). bond("leukotrieneA4",single,5,6). bond("leukotrieneA4",single,6,7).
bond("leukotrieneA4",double,7,8). bond("leukotrieneA4",single,8,9). bond("leukotrieneA4",double,9,10).
bond("leukotrieneA4",single,10,11). bond("leukotrieneA4",double,11,12). bond("leukotrieneA4",single,12,13).
bond("leukotrieneA4",single,13,14). bond("leukotrieneA4",double,14,15). bond("leukotrieneA4",single,15,16).
bond("leukotrieneA4",single,16,17). bond("leukotrieneA4",single,17,18). bond("leukotrieneA4",single,18,19).
bond("leukotrieneA4",single,19,20). bond("leukotrieneA4",single,1,21). bond("leukotrieneA4",double,1,22).
bond("leukotrieneA4",single,5,23). bond("leukotrieneA4",single,6,23).

% leukotrieneC4
atom("leukotrieneC4",1..20,carb). atom("leukotrieneC4",21..23,oxyg). atom("leukotrieneC4",24,glut).
bond("leukotrieneC4",single,1,2). bond("leukotrieneC4",single,2,3). bond("leukotrieneC4",single,3,4).
bond("leukotrieneC4",single,4,5). bond("leukotrieneC4",single,5,6). bond("leukotrieneC4",single,6,7).
bond("leukotrieneC4",double,7,8). bond("leukotrieneC4",single,8,9). bond("leukotrieneC4",double,9,10).
bond("leukotrieneC4",single,10,11). bond("leukotrieneC4",double,11,12). bond("leukotrieneC4",single,12,13).
bond("leukotrieneC4",single,13,14). bond("leukotrieneC4",double,14,15). bond("leukotrieneC4",single,15,16).
bond("leukotrieneC4",single,16,17). bond("leukotrieneC4",single,17,18). bond("leukotrieneC4",single,18,19).
bond("leukotrieneC4",single,19,20). bond("leukotrieneC4",single,1,21). bond("leukotrieneC4",double,1,22).
bond("leukotrieneC4",single,5,23). bond("leukotrieneC4",single,6,24).


%Shared domain.

atomDomain(fatty_acid,1..18,carb). atomDomain(fatty_acid,19..20,oxyg).
bondDomain(fatty_acid,single,1,2). bondDomain(fatty_acid,single,2,3). bondDomain(fatty_acid,single,3,4).
bondDomain(fatty_acid,single,4,5). bondDomain(fatty_acid,single,5,6). bondDomain(fatty_acid,single,6,7).
bondDomain(fatty_acid,single,7,8). bondDomain(fatty_acid,single,8,9). bondDomain(fatty_acid,single,9,10).
bondDomain(fatty_acid,single,10,11). bondDomain(fatty_acid,single,11,12). bondDomain(fatty_acid,single,12,13).
bondDomain(fatty_acid,single,13,14). bondDomain(fatty_acid,single,14,15). bondDomain(fatty_acid,single,15,16).
bondDomain(fatty_acid,single,16,17). bondDomain(fatty_acid,single,17,18). bondDomain(fatty_acid,single,1,19).
bondDomain(fatty_acid,double,1,20).

%* Known reactions
%%%%%%%%%%%%%%%%%%%
reaction(Reaction type (name), Compound 1 [, Active site atom number], [Compound 2 , Active site atom number], Resulting compound)
*%

% Reactions in Laminaria digitata.

reaction(lipoxygenase_rxn, "linoleic_acid", "13-HpODE").
reaction(epoxidation_L, "13-HpODE", "12,13-EpOME").

% Reactions in Ectocarpus sp7.

reaction(rxn_8495, "linoleic_acid", "9-HPOD").
reaction(epoxidation_E, "9-HPOD", "9,10-epoxy-11-hydroxy-12-octdecenoic").

% Reactions from animals

reaction(leukotriene-c4-synthase-rxn, "leukotrieneA4", "leukotrieneC4").

% Molecule Source.

source("linoleic_acid").


% Initiation and goal of incremental grounding.
goal(pathway("linoleic_acid","12-hydroxy-13-glutathion-OME")).

%Bionavi mol
%%%%%%%%%%%%

%b320
atom("b320",1,nitr).
atom("b320",7,nitr).
atom("b320",13,nitr).
atom("b320",2..5,carb).
atom("b320",8..9,carb).
atom("b320",11,carb).
atom("b320",14..15,carb).
atom("b320",18,carb).
atom("b320",6,oxyg).
atom("b320",12,oxyg).
atom("b320",16..17,oxyg).
atom("b320",19..20,oxyg).
atom("b320",10,sulf).
bond("b320",single,1,2).
bond("b320",single,2,3).
bond("b320",single,3,4).
bond("b320",single,4,5).
bond("b320",double,5,6).
bond("b320",single,5,7).
bond("b320",single,7,8).
bond("b320",single,8,9).
bond("b320",single,9,10).
bond("b320",single,8,11).
bond("b320",double,11,12).
bond("b320",single,11,13).
bond("b320",single,13,14).
bond("b320",single,14,15).
bond("b320",double,15,16).
bond("b320",single,15,17).
bond("b320",single,2,18).
bond("b320",double,18,19).
bond("b320",single,18,20).

%m0
atom("m0",1..6,carb).
atom("m0",8..9,carb).
atom("m0",11,carb).
atom("m0",13..15,carb).
atom("m0",17,carb).
atom("m0",20,carb).
atom("m0",23..24,carb).
atom("m0",27,carb).
atom("m0",29..39,carb).
atom("m0",7,sulf).
atom("m0",10,nitr).
atom("m0",16,nitr).
atom("m0",22,nitr).
atom("m0",12,oxyg).
atom("m0",18..19,oxyg).
atom("m0",21,oxyg).
atom("m0",25..26,oxyg).
atom("m0",28,oxyg).
atom("m0",40..41,oxyg).
bond("m0",single,1,2).
bond("m0",single,2,3).
bond("m0",single,3,4).
bond("m0",single,4,5).
bond("m0",single,5,6).
bond("m0",single,6,7).
bond("m0",single,7,8).
bond("m0",single,8,9).
bond("m0",single,9,10).
bond("m0",single,10,11).
bond("m0",double,11,12).
bond("m0",single,11,13).
bond("m0",single,13,14).
bond("m0",single,14,15).
bond("m0",single,15,16).
bond("m0",single,15,17).
bond("m0",double,17,18).
bond("m0",single,17,19).
bond("m0",single,9,20).
bond("m0",double,20,21).
bond("m0",single,20,22).
bond("m0",single,22,23).
bond("m0",single,23,24).
bond("m0",double,24,25).
bond("m0",single,24,26).
bond("m0",single,6,27).
bond("m0",single,27,28).
bond("m0",single,27,29).
bond("m0",single,29,30).
bond("m0",double,30,31).
bond("m0",single,31,32).
bond("m0",single,32,33).
bond("m0",single,33,34).
bond("m0",single,34,35).
bond("m0",single,35,36).
bond("m0",single,36,37).
bond("m0",single,37,38).
bond("m0",single,38,39).
bond("m0",double,39,40).
bond("m0",single,39,41).

%m2
atom("m2",1..6,carb).
atom("m2",8..19,carb).
atom("m2",7,oxyg).
atom("m2",20..21,oxyg).
bond("m2",single,1,2).
bond("m2",single,2,3).
bond("m2",single,3,4).
bond("m2",single,4,5).
bond("m2",single,5,6).
bond("m2",single,6,7).
bond("m2",single,7,8).
bond("m2",single,8,9).
bond("m2",single,9,10).
bond("m2",double,10,11).
bond("m2",single,11,12).
bond("m2",single,12,13).
bond("m2",single,13,14).
bond("m2",single,14,15).
bond("m2",single,15,16).
bond("m2",single,16,17).
bond("m2",single,17,18).
bond("m2",single,18,19).
bond("m2",double,19,20).
bond("m2",single,19,21).
bond("m2",single,8,6).

%b103721
atom("b103721",1..18,carb).
atom("b103721",19..20,oxyg).
bond("b103721",single,1,2).
bond("b103721",single,2,3).
bond("b103721",single,3,4).
bond("b103721",single,4,5).
bond("b103721",single,5,6).
bond("b103721",double,6,7).
bond("b103721",single,7,8).
bond("b103721",single,8,9).
bond("b103721",double,9,10).
bond("b103721",single,10,11).
bond("b103721",single,11,12).
bond("b103721",single,12,13).
bond("b103721",single,13,14).
bond("b103721",single,14,15).
bond("b103721",single,15,16).
bond("b103721",single,16,17).
bond("b103721",single,17,18).
bond("b103721",double,18,19).
bond("b103721",single,18,20).

%b103731
atom("b103731",1,nitr).
atom("b103731",2..3,carb).
atom("b103731",4..5,oxyg).
bond("b103731",single,1,2).
bond("b103731",single,2,3).
bond("b103731",double,3,4).
bond("b103731",single,3,5).

%m5
atom("m5",1..6,carb).
atom("m5",8..9,carb).
atom("m5",11,carb).
atom("m5",13..15,carb).
atom("m5",17,carb).
atom("m5",20,carb).
atom("m5",23,carb).
atom("m5",25..35,carb).
atom("m5",7,sulf).
atom("m5",10,nitr).
atom("m5",16,nitr).
atom("m5",12,oxyg).
atom("m5",18..19,oxyg).
atom("m5",21..22,oxyg).
atom("m5",24,oxyg).
atom("m5",36..37,oxyg).
bond("m5",single,1,2).
bond("m5",single,2,3).
bond("m5",single,3,4).
bond("m5",single,4,5).
bond("m5",single,5,6).
bond("m5",single,6,7).
bond("m5",single,7,8).
bond("m5",single,8,9).
bond("m5",single,9,10).
bond("m5",single,10,11).
bond("m5",double,11,12).
bond("m5",single,11,13).
bond("m5",single,13,14).
bond("m5",single,14,15).
bond("m5",single,15,16).
bond("m5",single,15,17).
bond("m5",double,17,18).
bond("m5",single,17,19).
bond("m5",single,9,20).
bond("m5",double,20,21).
bond("m5",single,20,22).
bond("m5",single,6,23).
bond("m5",single,23,24).
bond("m5",single,23,25).
bond("m5",single,25,26).
bond("m5",double,26,27).
bond("m5",single,27,28).
bond("m5",single,28,29).
bond("m5",single,29,30).
bond("m5",single,30,31).
bond("m5",single,31,32).
bond("m5",single,32,33).
bond("m5",single,33,34).
bond("m5",single,34,35).
bond("m5",double,35,36).
bond("m5",single,35,37).

%b227
atom("b227",1,nitr).
atom("b227",7,nitr).
atom("b227",2..5,carb).
atom("b227",8..9,carb).
atom("b227",11,carb).
atom("b227",14,carb).
atom("b227",6,oxyg).
atom("b227",12..13,oxyg).
atom("b227",15..16,oxyg).
atom("b227",10,sulf).
bond("b227",single,1,2).
bond("b227",single,2,3).
bond("b227",single,3,4).
bond("b227",single,4,5).
bond("b227",double,5,6).
bond("b227",single,5,7).
bond("b227",single,7,8).
bond("b227",single,8,9).
bond("b227",single,9,10).
bond("b227",single,8,11).
bond("b227",double,11,12).
bond("b227",single,11,13).
bond("b227",single,2,14).
bond("b227",double,14,15).
bond("b227",single,14,16).

%b2
atom("b2",1,nitr).
atom("b2",2..5,carb).
atom("b2",8,carb).
atom("b2",6..7,oxyg).
atom("b2",9..10,oxyg).
bond("b2",single,1,2).
bond("b2",single,2,3).
bond("b2",single,3,4).
bond("b2",single,4,5).
bond("b2",double,5,6).
bond("b2",single,5,7).
bond("b2",single,2,8).
bond("b2",double,8,9).
bond("b2",single,8,10).

%m8
atom("m8",1..6,carb).
atom("m8",8..9,carb).
atom("m8",11,carb).
atom("m8",14,carb).
atom("m8",16..26,carb).
atom("m8",7,sulf).
atom("m8",10,nitr).
atom("m8",12..13,oxyg).
atom("m8",15,oxyg).
atom("m8",27..28,oxyg).
bond("m8",single,1,2).
bond("m8",single,2,3).
bond("m8",single,3,4).
bond("m8",single,4,5).
bond("m8",single,5,6).
bond("m8",single,6,7).
bond("m8",single,7,8).
bond("m8",single,8,9).
bond("m8",single,9,10).
bond("m8",single,9,11).
bond("m8",double,11,12).
bond("m8",single,11,13).
bond("m8",single,6,14).
bond("m8",single,14,15).
bond("m8",single,14,16).
bond("m8",single,16,17).
bond("m8",double,17,18).
bond("m8",single,18,19).
bond("m8",single,19,20).
bond("m8",single,20,21).
bond("m8",single,21,22).
bond("m8",single,22,23).
bond("m8",single,23,24).
bond("m8",single,24,25).
bond("m8",single,25,26).
bond("m8",double,26,27).
bond("m8",single,26,28).

%b103732
atom("b103732",1,nitr).
atom("b103732",2..3,carb).
atom("b103732",5,carb).
atom("b103732",4,sulf).
atom("b103732",6..7,oxyg).
bond("b103732",single,1,2).
bond("b103732",single,2,3).
bond("b103732",single,3,4).
bond("b103732",single,2,5).
bond("b103732",double,5,6).
bond("b103732",single,5,7).

%m10
atom("m10",1..6,carb).
atom("m10",8,carb).
atom("m10",10..20,carb).
atom("m10",7,sulf).
atom("m10",9,oxyg).
atom("m10",21..22,oxyg).
bond("m10",single,1,2).
bond("m10",single,2,3).
bond("m10",single,3,4).
bond("m10",single,4,5).
bond("m10",single,5,6).
bond("m10",single,6,7).
bond("m10",single,6,8).
bond("m10",single,8,9).
bond("m10",single,8,10).
bond("m10",single,10,11).
bond("m10",double,11,12).
bond("m10",single,12,13).
bond("m10",single,13,14).
bond("m10",single,14,15).
bond("m10",single,15,16).
bond("m10",single,16,17).
bond("m10",single,17,18).
bond("m10",single,18,19).
bond("m10",single,19,20).
bond("m10",double,20,21).
bond("m10",single,20,22).

%m11
atom("m11",1..6,carb).
atom("m11",8..19,carb).
atom("m11",7,oxyg).
atom("m11",20..21,oxyg).
bond("m11",single,1,2).
bond("m11",single,2,3).
bond("m11",single,3,4).
bond("m11",single,4,5).
bond("m11",single,5,6).
bond("m11",single,6,7).
bond("m11",single,7,8).
bond("m11",single,8,9).
bond("m11",single,9,10).
bond("m11",double,10,11).
bond("m11",single,11,12).
bond("m11",single,12,13).
bond("m11",single,13,14).
bond("m11",single,14,15).
bond("m11",single,15,16).
bond("m11",single,16,17).
bond("m11",single,17,18).
bond("m11",single,18,19).
bond("m11",double,19,20).
bond("m11",single,19,21).
bond("m11",single,8,6).

%m12
atom("m12",1..6,carb).
atom("m12",8..9,carb).
atom("m12",11,carb).
atom("m12",13..15,carb).
atom("m12",17,carb).
atom("m12",20,carb).
atom("m12",23,carb).
atom("m12",25..36,carb).
atom("m12",7,sulf).
atom("m12",10,nitr).
atom("m12",16,nitr).
atom("m12",12,oxyg).
atom("m12",18..19,oxyg).
atom("m12",21..22,oxyg).
atom("m12",24,oxyg).
atom("m12",37..38,oxyg).
bond("m12",single,1,2).
bond("m12",single,2,3).
bond("m12",single,3,4).
bond("m12",single,4,5).
bond("m12",single,5,6).
bond("m12",single,6,7).
bond("m12",single,7,8).
bond("m12",single,8,9).
bond("m12",single,9,10).
bond("m12",single,10,11).
bond("m12",double,11,12).
bond("m12",single,11,13).
bond("m12",single,13,14).
bond("m12",single,14,15).
bond("m12",single,15,16).
bond("m12",single,15,17).
bond("m12",double,17,18).
bond("m12",single,17,19).
bond("m12",single,9,20).
bond("m12",double,20,21).
bond("m12",single,20,22).
bond("m12",single,6,23).
bond("m12",single,23,24).
bond("m12",single,23,25).
bond("m12",single,25,26).
bond("m12",double,26,27).
bond("m12",single,27,28).
bond("m12",single,28,29).
bond("m12",single,29,30).
bond("m12",single,30,31).
bond("m12",single,31,32).
bond("m12",single,32,33).
bond("m12",single,33,34).
bond("m12",single,34,35).
bond("m12",single,35,36).
bond("m12",double,36,37).
bond("m12",single,36,38).

%m13
atom("m13",1..6,carb).
atom("m13",8..9,carb).
atom("m13",11,carb).
atom("m13",13..15,carb).
atom("m13",17,carb).
atom("m13",20,carb).
atom("m13",23..24,carb).
atom("m13",27,carb).
atom("m13",29..39,carb).
atom("m13",7,sulf).
atom("m13",10,nitr).
atom("m13",16,nitr).
atom("m13",22,nitr).
atom("m13",12,oxyg).
atom("m13",18..19,oxyg).
atom("m13",21,oxyg).
atom("m13",25..26,oxyg).
atom("m13",28,oxyg).
atom("m13",40..41,oxyg).
bond("m13",single,1,2).
bond("m13",single,2,3).
bond("m13",single,3,4).
bond("m13",single,4,5).
bond("m13",single,5,6).
bond("m13",single,6,7).
bond("m13",single,7,8).
bond("m13",single,8,9).
bond("m13",single,9,10).
bond("m13",single,10,11).
bond("m13",double,11,12).
bond("m13",single,11,13).
bond("m13",single,13,14).
bond("m13",single,14,15).
bond("m13",single,15,16).
bond("m13",single,15,17).
bond("m13",double,17,18).
bond("m13",single,17,19).
bond("m13",single,9,20).
bond("m13",double,20,21).
bond("m13",single,20,22).
bond("m13",single,22,23).
bond("m13",single,23,24).
bond("m13",double,24,25).
bond("m13",single,24,26).
bond("m13",single,6,27).
bond("m13",double,27,28).
bond("m13",single,27,29).
bond("m13",single,29,30).
bond("m13",double,30,31).
bond("m13",single,31,32).
bond("m13",single,32,33).
bond("m13",single,33,34).
bond("m13",single,34,35).
bond("m13",single,35,36).
bond("m13",single,36,37).
bond("m13",single,37,38).
bond("m13",single,38,39).
bond("m13",double,39,40).
bond("m13",single,39,41).

%m14
atom("m14",1..6,carb).
atom("m14",8..19,carb).
atom("m14",7,oxyg).
atom("m14",20..21,oxyg).
bond("m14",single,1,2).
bond("m14",single,2,3).
bond("m14",single,3,4).
bond("m14",single,4,5).
bond("m14",single,5,6).
bond("m14",single,6,7).
bond("m14",single,7,8).
bond("m14",single,8,9).
bond("m14",single,9,10).
bond("m14",double,10,11).
bond("m14",single,11,12).
bond("m14",single,12,13).
bond("m14",single,13,14).
bond("m14",single,14,15).
bond("m14",single,15,16).
bond("m14",single,16,17).
bond("m14",single,17,18).
bond("m14",single,18,19).
bond("m14",double,19,20).
bond("m14",single,19,21).
bond("m14",single,8,6).

%Bionavi rxn
%%%%%%%%%%%%
