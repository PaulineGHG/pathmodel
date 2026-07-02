%*
REACTIONS
====================================================================================================
*%

% Reaction1
reactant("Reaction1","Reaction1_reactant1","False").
reactant("Reaction1","Reaction1_reactant2","False").
product("Reaction1","Reaction1_product1","False").

% Reaction2
reactant("Reaction2","Reaction2_reactant1","False").
product("Reaction2","Reaction2_product1","False").

%*
CHEMICALS
====================================================================================================
*%

% Source
atom("Source",1,c).
atom("Source",2,c).
atom("Source",3,c).
atom("Source",4,o).
atom("Source",5,c).
atom("Source",6,c).
atom("Source",7,c).
atom("Source",8,c).
atom("Source",9,c).
atom("Source",10,o).
atom("Source",11,o).
bond("Source",double,1,2).
bond("Source",single,2,3).
bond("Source",single,3,4).
bond("Source",single,4,5).
bond("Source",single,5,6).
bond("Source",single,6,7).
bond("Source",single,7,8).
bond("Source",single,8,9).
bond("Source",double,9,10).
bond("Source",single,9,11).
bond("Source",single,5,3).

% Target
atom("Target",1,c).
atom("Target",2,c).
atom("Target",3,c).
atom("Target",4,o).
atom("Target",5,c).
atom("Target",6,c).
atom("Target",7,c).
atom("Target",8,c).
atom("Target",9,c).
atom("Target",10,o).
atom("Target",11,o).
atom("Target",12,s).
atom("Target",13,c).
atom("Target",14,o).
atom("Target",15,o).
bond("Target",double,1,2).
bond("Target",single,2,3).
bond("Target",single,3,4).
bond("Target",single,3,5).
bond("Target",single,5,6).
bond("Target",single,6,7).
bond("Target",single,7,8).
bond("Target",single,8,9).
bond("Target",double,9,10).
bond("Target",single,9,11).
bond("Target",single,5,12).
bond("Target",single,12,13).
bond("Target",double,13,14).
bond("Target",single,13,15).

% Reaction1_reactant1
atom("Reaction1_reactant1",1,c).
atom("Reaction1_reactant1",2,c).
atom("Reaction1_reactant1",3,c).
atom("Reaction1_reactant1",4,c).
atom("Reaction1_reactant1",5,o).
atom("Reaction1_reactant1",6,c).
atom("Reaction1_reactant1",7,c).
atom("Reaction1_reactant1",8,c).
atom("Reaction1_reactant1",9,o).
atom("Reaction1_reactant1",10,o).
bond("Reaction1_reactant1",single,1,2).
bond("Reaction1_reactant1",double,2,3).
bond("Reaction1_reactant1",single,3,4).
bond("Reaction1_reactant1",single,4,5).
bond("Reaction1_reactant1",single,5,6).
bond("Reaction1_reactant1",single,6,7).
bond("Reaction1_reactant1",single,7,8).
bond("Reaction1_reactant1",double,8,9).
bond("Reaction1_reactant1",single,8,10).
bond("Reaction1_reactant1",single,6,4).

% Reaction1_reactant2
atom("Reaction1_reactant2",1,s).
atom("Reaction1_reactant2",2,c).
atom("Reaction1_reactant2",3,o).
atom("Reaction1_reactant2",4,o).
bond("Reaction1_reactant2",single,1,2).
bond("Reaction1_reactant2",double,2,3).
bond("Reaction1_reactant2",single,2,4).

% Reaction1_product1
atom("Reaction1_product1",1,c).
atom("Reaction1_product1",2,c).
atom("Reaction1_product1",3,c).
atom("Reaction1_product1",4,c).
atom("Reaction1_product1",5,o).
atom("Reaction1_product1",6,c).
atom("Reaction1_product1",7,c).
atom("Reaction1_product1",8,c).
atom("Reaction1_product1",9,o).
atom("Reaction1_product1",10,o).
atom("Reaction1_product1",11,s).
atom("Reaction1_product1",12,c).
atom("Reaction1_product1",13,o).
atom("Reaction1_product1",14,o).
bond("Reaction1_product1",single,1,2).
bond("Reaction1_product1",double,2,3).
bond("Reaction1_product1",single,3,4).
bond("Reaction1_product1",single,4,5).
bond("Reaction1_product1",single,4,6).
bond("Reaction1_product1",single,6,7).
bond("Reaction1_product1",single,7,8).
bond("Reaction1_product1",double,8,9).
bond("Reaction1_product1",single,8,10).
bond("Reaction1_product1",single,6,11).
bond("Reaction1_product1",single,11,12).
bond("Reaction1_product1",double,12,13).
bond("Reaction1_product1",single,12,14).

% Reaction2_reactant1
atom("Reaction2_reactant1",1,c).
atom("Reaction2_reactant1",2,c).
atom("Reaction2_reactant1",3,c).
atom("Reaction2_reactant1",4,o).
atom("Reaction2_reactant1",5,o).
atom("Reaction2_reactant1",6,c).
atom("Reaction2_reactant1",7,c).
atom("Reaction2_reactant1",8,c).
atom("Reaction2_reactant1",9,c).
atom("Reaction2_reactant1",10,o).
atom("Reaction2_reactant1",11,o).
bond("Reaction2_reactant1",single,1,2).
bond("Reaction2_reactant1",single,2,3).
bond("Reaction2_reactant1",single,3,4).
bond("Reaction2_reactant1",single,4,5).
bond("Reaction2_reactant1",single,3,6).
bond("Reaction2_reactant1",double,6,7).
bond("Reaction2_reactant1",single,7,8).
bond("Reaction2_reactant1",single,8,9).
bond("Reaction2_reactant1",double,9,10).
bond("Reaction2_reactant1",single,9,11).

% Reaction2_product1
atom("Reaction2_product1",1,c).
atom("Reaction2_product1",2,c).
atom("Reaction2_product1",3,c).
atom("Reaction2_product1",4,o).
atom("Reaction2_product1",5,o).
atom("Reaction2_product1",6,c).
atom("Reaction2_product1",7,c).
atom("Reaction2_product1",8,c).
atom("Reaction2_product1",9,c).
atom("Reaction2_product1",10,o).
atom("Reaction2_product1",11,o).
bond("Reaction2_product1",single,1,2).
bond("Reaction2_product1",single,2,3).
bond("Reaction2_product1",single,3,4).
bond("Reaction2_product1",single,3,5).
bond("Reaction2_product1",single,5,6).
bond("Reaction2_product1",double,6,7).
bond("Reaction2_product1",single,7,8).
bond("Reaction2_product1",single,8,9).
bond("Reaction2_product1",double,9,10).
bond("Reaction2_product1",single,9,11).

%*
SOURCE - GOAL
====================================================================================================
*%

% SOURCE
source("Source").

% GOAL
goal(pathway("Source","Target")).

%*
MAPPINGS
====================================================================================================
*%

% Reaction1
atomMappingReactant("Reaction1",c,1).
atomMappingReactant("Reaction1",c,2).
atomMappingReactant("Reaction1",c,3).
atomMappingReactant("Reaction1",c,4).
atomMappingReactant("Reaction1",o,5).
atomMappingReactant("Reaction1",c,6).
atomMappingReactant("Reaction1",c,7).
atomMappingReactant("Reaction1",c,8).
atomMappingReactant("Reaction1",o,9).
atomMappingReactant("Reaction1",o,10).
bondMappingReactant("Reaction1",1,2,single).
bondMappingReactant("Reaction1",2,3,double).
bondMappingReactant("Reaction1",3,4,single).
bondMappingReactant("Reaction1",4,5,single).
bondMappingReactant("Reaction1",5,6,single).
bondMappingReactant("Reaction1",6,7,single).
bondMappingReactant("Reaction1",7,8,single).
bondMappingReactant("Reaction1",8,9,double).
bondMappingReactant("Reaction1",8,10,single).
bondMappingReactant("Reaction1",6,4,single).
atomMappingReactant("Reaction1",s,11).
atomMappingReactant("Reaction1",c,12).
atomMappingReactant("Reaction1",o,13).
atomMappingReactant("Reaction1",o,14).
bondMappingReactant("Reaction1",11,12,single).
bondMappingReactant("Reaction1",12,13,double).
bondMappingReactant("Reaction1",12,14,single).
atomMappingProduct("Reaction1",c,1).
atomMappingProduct("Reaction1",c,2).
atomMappingProduct("Reaction1",c,3).
atomMappingProduct("Reaction1",c,4).
atomMappingProduct("Reaction1",o,5).
atomMappingProduct("Reaction1",c,6).
atomMappingProduct("Reaction1",c,7).
atomMappingProduct("Reaction1",c,8).
atomMappingProduct("Reaction1",o,9).
atomMappingProduct("Reaction1",o,10).
atomMappingProduct("Reaction1",s,11).
atomMappingProduct("Reaction1",c,12).
atomMappingProduct("Reaction1",o,13).
atomMappingProduct("Reaction1",o,14).
bondMappingProduct("Reaction1",1,2,single).
bondMappingProduct("Reaction1",2,3,double).
bondMappingProduct("Reaction1",3,4,single).
bondMappingProduct("Reaction1",4,5,single).
bondMappingProduct("Reaction1",4,6,single).
bondMappingProduct("Reaction1",6,7,single).
bondMappingProduct("Reaction1",7,8,single).
bondMappingProduct("Reaction1",8,9,double).
bondMappingProduct("Reaction1",8,10,single).
bondMappingProduct("Reaction1",6,11,single).
bondMappingProduct("Reaction1",11,12,single).
bondMappingProduct("Reaction1",12,13,double).
bondMappingProduct("Reaction1",12,14,single).

% Reaction2
atomMappingReactant("Reaction2",c,1).
atomMappingReactant("Reaction2",c,2).
atomMappingReactant("Reaction2",c,3).
atomMappingReactant("Reaction2",o,5).
atomMappingReactant("Reaction2",o,4).
atomMappingReactant("Reaction2",c,6).
atomMappingReactant("Reaction2",c,7).
atomMappingReactant("Reaction2",c,8).
atomMappingReactant("Reaction2",c,9).
atomMappingReactant("Reaction2",o,10).
atomMappingReactant("Reaction2",o,11).
bondMappingReactant("Reaction2",1,2,single).
bondMappingReactant("Reaction2",2,3,single).
bondMappingReactant("Reaction2",3,5,single).
bondMappingReactant("Reaction2",5,4,single).
bondMappingReactant("Reaction2",3,6,single).
bondMappingReactant("Reaction2",6,7,double).
bondMappingReactant("Reaction2",7,8,single).
bondMappingReactant("Reaction2",8,9,single).
bondMappingReactant("Reaction2",9,10,double).
bondMappingReactant("Reaction2",9,11,single).
atomMappingProduct("Reaction2",c,1).
atomMappingProduct("Reaction2",c,2).
atomMappingProduct("Reaction2",c,3).
atomMappingProduct("Reaction2",o,4).
atomMappingProduct("Reaction2",o,5).
atomMappingProduct("Reaction2",c,6).
atomMappingProduct("Reaction2",c,7).
atomMappingProduct("Reaction2",c,8).
atomMappingProduct("Reaction2",c,9).
atomMappingProduct("Reaction2",o,10).
atomMappingProduct("Reaction2",o,11).
bondMappingProduct("Reaction2",1,2,single).
bondMappingProduct("Reaction2",2,3,single).
bondMappingProduct("Reaction2",3,4,single).
bondMappingProduct("Reaction2",3,5,single).
bondMappingProduct("Reaction2",5,6,single).
bondMappingProduct("Reaction2",6,7,double).
bondMappingProduct("Reaction2",7,8,single).
bondMappingProduct("Reaction2",8,9,single).
bondMappingProduct("Reaction2",9,10,double).
bondMappingProduct("Reaction2",9,11,single).
