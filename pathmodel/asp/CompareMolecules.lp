 %*
Search difference between molecules occurring in the same reaction.
This script is used to check if reactions are well encoded.
*%

%* Test to see difference between two molecules
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*%

% Pour chaque reaction d'un reactant et d'un produit on crée une comparaison de ce reactant et du produit pour cette reaction
compareMolecules(ReactionName, Reactant, Product):- reaction(ReactionName, Reactant, Product).

% Pour chaque comparaison entre un réactant et un produit on crée une comparaison entre le produit et le réactant pour la même réaction
compareMolecules(ReactionName,Reactant, Product):- compareMolecules(ReactionName, Product, Reactant).

%* 
Pour chaque comparaison entre 2 molécules, si il existe des liaisons entre deux même numéro d'atome non présents chez les deux molécules, on crée une différence
de ces deux molécules pour cette liaison, ces atomes, ces molécules et la reaction associée --> Equivaut à une transformation
*%
diffMolecules(ReactionName, MoleculeName1, MoleculeName2, BondType, FirstAtomNumber, SecondAtomNumber):-
		    bond(MoleculeName1, BondType, FirstAtomNumber, SecondAtomNumber),
		not bond(MoleculeName2, BondType, FirstAtomNumber, SecondAtomNumber), compareMolecules(ReactionName, MoleculeName1, MoleculeName2).


% Affiche toutes les différences (=éléments d'une transformation) (DiffMolecules ayant 6 arguments)
#show diffMolecules/6.