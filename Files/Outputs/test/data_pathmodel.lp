% Molecule 1
atom("molecule_1",1..4,carb). 
bond("molecule_1",single,1,2). bond("molecule_1",single,1,3). bond("molecule_1",single,2,3).
bond("molecule_1",single,2,4).

% Molecule 2
atom("molecule_2",1..4,carb). 
bond("molecule_2",single,1,2). bond("molecule_2",single,1,3). bond("molecule_2",single,2,3).
bond("molecule_2",double,2,4).

% Molecule 3
atom("molecule_3",1..6,carb). 
bond("molecule_3",single,1,2). bond("molecule_3",single,1,3). bond("molecule_3",single,1,6).
bond("molecule_3",single,2,3). bond("molecule_3",single,2,4). bond("molecule_3",single,3,6).
bond("molecule_3",single,5,6).

% Molecule 4
atom("molecule_4",1..6,carb). 
bond("molecule_4",single,1,2). bond("molecule_4",single,1,3). bond("molecule_4",single,1,6).
bond("molecule_4",single,2,3). bond("molecule_4",double,2,4). bond("molecule_4",single,3,6).
bond("molecule_4",single,5,6).

% Molecule 5
atom("molecule_5",1..7,carb). 
bond("molecule_5",single,1,2). bond("molecule_5",single,1,3). bond("molecule_5",single,1,6).
bond("molecule_5",single,1,7). bond("molecule_5",single,2,3). bond("molecule_5",single,2,4).
bond("molecule_5",double,3,6). bond("molecule_5",single,5,6).

% Shared domain.

atomDomain(triangle,1..3,carb).
bondDomain(triangle,single,1,2). bondDomain(triangle,single,1,3).
bondDomain(triangle,single,2,3).

% Known reaction.

reaction(reduction, "molecule_1", "molecule_2").

% Molecule Source.

source("molecule_3").
source("molecule_5").

% Known MZ.

mzfiltering(921341).

% Initiation and goal of incremental grounding.
init(pathway("molecule_1","molecule_2")).
goal(pathway("molecule_3","molecule_4")).

step_limit(100).
domain("molecule_1",triangle).
domain("molecule_2",triangle).
domain("molecule_3",triangle).
domain("molecule_4",triangle).
domain("molecule_5",triangle).
moleculeComposition("molecule_1",4,8,0,0,0).
moleculeComposition("molecule_2",4,6,0,0,0).
moleculeComposition("molecule_3",6,10,0,0,0).
moleculeComposition("molecule_4",6,8,0,0,0).
moleculeComposition("molecule_5",7,10,0,0,0).
moleculeMZ("molecule_1",561020).
moleculeMZ("molecule_2",540872).
moleculeMZ("molecule_3",821382).
moleculeMZ("molecule_4",801234).
moleculeMZ("molecule_5",941489).
moleculeNbAtoms("molecule_1",12).
moleculeNbAtoms("molecule_1",4).
moleculeNbAtoms("molecule_2",10).
moleculeNbAtoms("molecule_2",4).
moleculeNbAtoms("molecule_3",16).
moleculeNbAtoms("molecule_3",6).
moleculeNbAtoms("molecule_4",14).
moleculeNbAtoms("molecule_4",6).
moleculeNbAtoms("molecule_5",17).
moleculeNbAtoms("molecule_5",7).
numberTotalBonds("molecule_1",4).
numberTotalBonds("molecule_2",4).
numberTotalBonds("molecule_3",7).
numberTotalBonds("molecule_4",7).
numberTotalBonds("molecule_5",8).
reactionMZ(reduction,20148).
diffBondAfterReaction(reduction,double,2,4).
diffBondBeforeReaction(reduction,single,2,4).
siteAfterReaction(reduction,"molecule_2").
siteAfterReaction(reduction,"molecule_4").
siteBeforeReaction(reduction,"molecule_1").
siteBeforeReaction(reduction,"molecule_3").
siteBeforeReaction(reduction,"molecule_5").
