%*
This script will create the reference transformation extracting a reference reaction center.
CHANGE ? -> To do this it compares links and atoms between two molecules that are implied in a reaction.
When it finds difference, the script will extract the links and atoms and write the result in a file.
*%

%*
Chemical variable:
valence(Atom type, Number of valence)
*%
valence(c, 4). valence(n, 3). valence(o, 2). valence(p, 5).

bond(Molecule, Atom2Position, Atom1Position, BondOrder) :- bond(Molecule, Atom1Position, Atom2Position, BondOrder).
bondMappingProduct(Molecule, Atom2Position, Atom1Position, BondOrder) :- bondMappingProduct(Molecule, Atom1Position, Atom2Position, BondOrder).
bondMappingReactant(Molecule, Atom2Position, Atom1Position, BondOrder) :- bondMappingReactant(Molecule, Atom1Position, Atom2Position, BondOrder).


%*
2. REACTANT / PRODUCT DIFFERENCES EXTRACTION
========================================================================================================================

The idea is to model a Condensed Graph of Reaction (CGR) extracting only the differences between the reactant and the
product as the reaction center :
What was added : atomAddition / bondAddition
What was removed : atomRemoval / bondRemoval
*%


    % 2.1 Compare atoms between molecules.
    % ------------------------------------

        % 2.1.1 Atoms added in the transformation = present in Product but not in Reactant.
atomDifference(Reaction, AtomElement, AtomPosition, addition) :- not atomMappingReactant(Reaction, AtomElement, AtomPosition),
                                                                 atomMappingProduct(Reaction, AtomElement, AtomPosition).

        % 2.1.2 Atoms removed in the transformation = present in Reactant but not in Product.
atomDifference(Reaction, AtomElement, AtomPosition, removal) :- atomMappingReactant(Reaction, AtomElement, AtomPosition),
                                                                not atomMappingProduct(Reaction, AtomElement, AtomPosition).


    % 2.2 Compare bonds between molecules.
    % ------------------------------------

        % 2.2.1 Bonds added in the transformation = present in Product but not in Reactant.
bondDifference(Reaction, Atom1Position, Atom2Position, BondOrder, addition) :- bondMappingProduct(Reaction, Atom1Position, Atom2Position, BondOrder),
                                                                               not bondMappingReactant(Reaction, Atom1Position, Atom2Position, BondOrder).

        % 2.2.2 Bonds removed in the transformation = present in Reactant but not in Product.
bondDifference(Reaction, Atom1Position, Atom2Position, BondOrder, removal) :- not bondMappingProduct(Reaction, Atom1Position, Atom2Position, BondOrder),
                                                                              bondMappingReactant(Reaction, Atom1Position, Atom2Position, BondOrder).


%*
3. Identification du bloc lié à la transformation
========================================================================================================================
L'objectif est d'intentifier le bloc du CGR regroupant les atomes et liaisons directement liées au centre de réaction
*%

%   3.1 Identifier les atomes du reactant impliqués dans la transformation
%   ----------------------------------------------------------------------

    % 3.1.1 REACTANT

        % Lien modifié entre ATOME et Atom2 dans le reactant
reactantAtomTransformed(Transformation, AtomElement, AtomPosition) :- bondDifference(Transformation, AtomPosition, Atom2Position, BondOrder, Modification),
                                                                      atomMappingReactant(Transformation, AtomElement, AtomPosition).


    % 3.1.2 PRODUCT

        % Lien modifié entre ATOME et Atom2 dans le produit
productAtomTransformed(Transformation, AtomElement, AtomPosition) :- bondDifference(Transformation, AtomPosition, Atom2Position, BondOrder, Modification),
                                                                     atomMappingProduct(Transformation, AtomElement, AtomPosition).



%   3.2 Identifier les liens liés à ces atomes
%   ------------------------------------------

    % 3.2.1 REACTANT

        % Lien entre un atome identifié et un autre atome dans le réactant
reactantTransformationCenterBond(Transformation, Atom1Position, Atom2Position, BondOrder) :- reactantAtomTransformed(Transformation, AtomElement, Atom1Position),
                                                                                             bondMappingReactant(Transformation, Atom1Position, Atom2Position, BondOrder).
        % Symétrie
reactantTransformationCenterBond(Transformation, Atom2Position, Atom1Position, BondOrder) :- reactantTransformationCenterBond(Transformation, Atom1Position, Atom2Position, BondOrder).

    % 3.2.2 PRODUCT

        % Lien entre un atome identifié et un autre atome dans le produit
productTransformationCenterBond(Transformation, Atom1Position, Atom2Position, BondOrder) :- productAtomTransformed(Transformation, AtomElement, Atom1Position),
                                                                                            bondMappingProduct(Transformation, Atom1Position, Atom2Position, BondOrder).
        % Symétrie
productTransformationCenterBond(Transformation, Atom2Position, Atom1Position, BondOrder) :- productTransformationCenterBond(Transformation, Atom1Position, Atom2Position, BondOrder).



%   3.3 Identifier les atomes liés à ces liens
%   ------------------------------------------

    % 3.3.1 REACTANT

reactantTransformationCenterAtom(Transformation, AtomElement, AtomPosition) :- reactantTransformationCenterBond(Transformation, AtomPosition, Atom2Position, BondOrder),
                                                                               atomMappingReactant(Transformation, AtomElement, AtomPosition).

    % 3.3.2 PRODUCT

productTransformationCenterAtom(Transformation, AtomElement, AtomPosition) :- productTransformationCenterBond(Transformation, AtomPosition, Atom2Position, BondOrder),
                                                                              atomMappingProduct(Transformation, AtomElement, AtomPosition).



#show atomDifference/4.
#show bondDifference/5.

#show reactantAtomTransformed/3.
#show productAtomTransformed/3.

#show reactantTransformationCenterBond/4.
#show reactantTransformationCenterAtom/3.

#show productTransformationCenterBond/4.
#show productTransformationCenterAtom/3.