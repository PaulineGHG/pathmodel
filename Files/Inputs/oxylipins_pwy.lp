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

% C18 fatty acid
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
