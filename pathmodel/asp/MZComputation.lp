%*
MZ Computation script
*%

%*
Chemical variable:
valence(Atom type, Number of valence)
*%
valence(carb, 4). valence(nitr, 3). valence(oxyg, 2). valence(phos, 5).

% Symmetric bonds.

%*
Si il y a une liaison pour une Molecule (A) donnée et un type de liaison (B) donné entre un 1er atome (C) et un second atome (D) alors on crée une
liaison 2 entre C et D.
Et si il y a une liaison entre D et C alors on crée une liaison 2 entre C et D aussi.
*%
bond2(A,B,C,D):- bond(A,B,C,D).
bond2(A,B,C,D):- bond(A,B,D,C).

% Si un input est déclaré pour une molécule (désignation de atom(Mol,_,_)), on en crée un objet metabolite à partir de son nom
metabolite(MoleculeName):- atom(MoleculeName,_,_).

%* Definition of atomic masses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*%

%*
Masses des principaux atomes
m_hydr = 1.0074, m_carb = 12.0107,  m_oxyg = 15.9994 ,  m_nitr=14.0067, m_phos = 30,9738.

Compute bonds between each atoms C,N,O, and add as many H as nb_H = valence (atom) - nb_bond(atom), with valence(oxyg)=2, valence(nitr)=3, valence(carb)=4 and valence(phos)=5(?).

mz(newcompound)=nb_carb(newcompound) x m_carb + nb_oxyg(newcompound) x m_oxyg + nb_nitr(newcompound) x m_nitr + nb_hydr(newcompound) x m_hydr + nb_phos(newcompound) x m_phos
*%

% M/Z Ratio
																
% numberHydrogens provides the number of hydrogens associated with each atoms.
% numberHydrogens(MoleculeName, AtomNumber, number of Hydrogen bonded with the atom)
numberHydrogens(MoleculeName, AtomNumber, NumberHydrogen) :-
	NumberOfBonds2={bond2(MoleculeName, double, AtomNumber, SecondAtomNumber)};
	NumberOfBonds1=#sum{1, SecondAtomNumber, BoundType: bond2(MoleculeName, BoundType, AtomNumber, SecondAtomNumber), BoundType != double};
	atom(MoleculeName, AtomNumber, AtomeType);
	valence(AtomeType, ValenceNumber);
	NumberHydrogen=ValenceNumber - NumberOfBonds1- 2*NumberOfBonds2; NumberHydrogen=0..ValenceNumber.

% moleculeComposition shows the number of Carbon, Hydrogen, Oxygen and Nitrogen in the molecule.
% moleculeComposition(MoleculeName, Number of Carbon, Number of Hydrogen, Number of Oxygen, Number of Nitrogen)
moleculeComposition(MoleculeName, NumberCarbon, NumberHydrogen, NumberOxygen, NumberNitrogen, NumberPhosphorus) :-
	NumberHydrogen=#sum{NumberHydrogenAtom, AtomNumber: numberHydrogens(MoleculeName, AtomNumber, NumberHydrogenAtom)},
	NumberCarbon={atom(MoleculeName, AtomNumber , carb)},
	NumberOxygen={atom(MoleculeName, AtomNumber , oxyg)},
	NumberNitrogen={atom(MoleculeName, AtomNumber , nitr)},
	NumberPhosphorus={atom(MoleculeName, AtomNumber , phos)},
	metabolite(MoleculeName).

% moleculeNbAtoms provides the total number of atoms in a compound.
moleculeNbAtoms(MoleculeName, NumberCarbon+NumberHydrogen+NumberOxygen+NumberNitrogen+NumberPhosphorus):-
	moleculeComposition(MoleculeName, NumberCarbon, NumberHydrogen, NumberOxygen, NumberNitrogen, NumberPhosphorus).

% moleculeNbAtoms provides the total number of atoms in a compound.
moleculeNbAtoms(MoleculeName, NumberCarbon + NumberOxygen + NumberNitrogen + NumberPhosphorus):-
	moleculeComposition(MoleculeName, NumberCarbon, NumberHydrogen, NumberOxygen, NumberNitrogen, NumberPhosphorus).

% numberTotalBonds provides the number of bond in a molecule.
numberTotalBonds(MoleculeName,NumberTotalOfBonds):-  NumberTotalOfBonds={bond(MoleculeName, BoundType, AtomNumber, SecondAtomNumber)}; metabolite(MoleculeName).

% moleculeMZ computes the M/Z ratio for each molecule.
% moleculeMZ(MoleculeName, M/Z ratio*10000)
% Because in ASP there is no decimal, all the atomic masses have been multiplied by 10 0000. If you want the real M/Z ratio divide the M/Z ratio by 10 0000.
% Approximation like at PubChem: 120110*NumberCarbon + 10080*NumberHydrogen + 159990*NumberOxygen + 140070*NumberNitrogen + 309738*NumberPhosphorus)
% 120107*NumberCarbon + 100794*NumberHydrogen + 159994*NumberOxygen + 140067*NumberNitrogen
% 1200960*NumberCarbon + 100784*NumberHydrogen + 1599903*NumberOxygen + 1400643*NumberNitrogen (here multiplied by 100 000) <- from https://www.degruyter.com/downloadpdf/j/pac.2016.88.issue-3/pac-2015-0305/pac-2015-0305.pdf
% 1201160*NumberCarbon + 100811*NumberHydrogen + 1599977*NumberOxygen + 1400728*NumberNitrogen (here multiplied by 100 000) <- from https://www.degruyter.com/downloadpdf/j/pac.2016.88.issue-3/pac-2015-0305/pac-2015-0305.pdf
moleculeMZ(MoleculeName, 120107*NumberCarbon + 10074*NumberHydrogen + 159994*NumberOxygen + 140067*NumberNitrogen + 309738*NumberPhosphorus):-
	moleculeComposition(MoleculeName, NumberCarbon, NumberHydrogen, NumberOxygen, NumberNitrogen, NumberPhosphorus).

% Compute MZ for reaction site.
reactionMZ(Reaction, ReactionMZ):- reaction(Reaction,Reactant,Product), moleculeMZ(Reactant,ReactantMZ), moleculeMZ(Product,ProductMZ),
									ProductMZ > ReactantMZ, ReactionMZ = ProductMZ - ReactantMZ.

reactionMZ(Reaction, ReactionMZ):- reaction(Reaction,Reactant,Product), moleculeMZ(Reactant,ReactantMZ), moleculeMZ(Product,ProductMZ),
									ProductMZ < ReactantMZ, ReactionMZ = ReactantMZ - ProductMZ.

%*
Definition of domain
domain(MoleculeName, DomainName)
*%

domain(MoleculeName, Domain) :- atomDomain(Domain,Atom,AtomType), bondDomain(Domain,_,Atom1,Atom2),
								atom(MoleculeName,Atom,AtomType), bond(MoleculeName,_,Atom1,Atom2).

#show domain/2.
#show moleculeComposition/6.
#show moleculeNbAtoms/2.
#show numberTotalBonds/2.
#show moleculeMZ/2.
#show reactionMZ/2.
