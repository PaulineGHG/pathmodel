%*
REACTIONS
====================================================================================================
*%

% Reaction1
reactant("Reaction1","M1","uni").
reactant("Reaction1","M2","uni").
product("Reaction1","M3","uni").


%*
CHEMICALS
====================================================================================================
*%

% Source
atom("Source",c,1).
atom("Source",c,2).
atom("Source",c,3).
atom("Source",c,4).
atom("Source",c,5).
atom("Source",c,6).
atom("Source",c,7).
atom("Source",c,8).
atom("Source",c,9).
atom("Source",o,10).
atom("Source",o,11).
bond("Source",1,2,single).
bond("Source",1,10,single).
bond("Source",1,9,double).
bond("Source",2,3,single).
bond("Source",3,4,single).
bond("Source",4,5,single).
bond("Source",5,6,single).
bond("Source",5,11,single).
bond("Source",6,11,single).
bond("Source",6,7,single).
bond("Source",7,8,double).

% Target
atom("Target",c,1).
atom("Target",c,2).
atom("Target",c,3).
atom("Target",c,4).
atom("Target",c,5).
atom("Target",c,6).
atom("Target",c,7).
atom("Target",c,8).
atom("Target",c,9).
atom("Target",s,10).
atom("Target",o,11).
atom("Target",o,12).
atom("Target",o,13).
atom("Target",o,14).
atom("Target",o,15).
bond("Target",1,2,single).
bond("Target",1,12,single).
bond("Target",1,11,double).
bond("Target",2,3,single).
bond("Target",3,4,single).
bond("Target",4,5,single).
bond("Target",5,10,single).
bond("Target",5,6,single).
bond("Target",6,13,single).
bond("Target",6,7,single).
bond("Target",7,8,double).
bond("Target",10,9,single).
bond("Target",9,15,single).
bond("Target",9,14,double).

% M1
atom("M1",c,1).
atom("M1",c,2).
atom("M1",c,3).
atom("M1",c,4).
atom("M1",c,5).
atom("M1",c,6).
atom("M1",c,7).
atom("M1",o,8).
atom("M1",o,9).
atom("M1",o,10).
bond("M1",1,2,single).
bond("M1",1,9,single).
bond("M1",1,8,double).
bond("M1",2,3,single).
bond("M1",3,4,single).
bond("M1",3,10,single).
bond("M1",4,10,single).
bond("M1",4,5,single).
bond("M1",5,6,double).
bond("M1",6,7,single).

% M2
atom("M2",c,1).
atom("M2",s,2).
atom("M2",o,3).
atom("M2",o,4).
bond("M2",1,2,single).
bond("M2",1,4,single).
bond("M2",1,3,double).

% M3
atom("M3",c,1).
atom("M3",c,2).
atom("M3",c,3).
atom("M3",c,4).
atom("M3",c,5).
atom("M3",c,6).
atom("M3",c,7).
atom("M3",c,8).
atom("M3",s,9).
atom("M3",o,10).
atom("M3",o,11).
atom("M3",o,12).
atom("M3",o,13).
atom("M3",o,14).
bond("M3",1,2,single).
bond("M3",1,13,single).
bond("M3",1,12,double).
bond("M3",2,3,single).
bond("M3",3,4,single).
bond("M3",4,14,single).
bond("M3",4,5,single).
bond("M3",5,6,double).
bond("M3",6,7,single).
bond("M3",3,9,single).
bond("M3",9,8,single).
bond("M3",8,11,single).
bond("M3",8,10,double).


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
atomMappingReactant("Reaction1",c,5).
atomMappingReactant("Reaction1",c,6).
atomMappingReactant("Reaction1",c,7).
atomMappingReactant("Reaction1",o,8).
atomMappingReactant("Reaction1",o,9).
atomMappingReactant("Reaction1",o,10).
bondMappingReactant("Reaction1",1,2,single).
bondMappingReactant("Reaction1",1,9,single).
bondMappingReactant("Reaction1",1,8,double).
bondMappingReactant("Reaction1",2,3,single).
bondMappingReactant("Reaction1",3,4,single).
bondMappingReactant("Reaction1",3,10,single).
bondMappingReactant("Reaction1",4,10,single).
bondMappingReactant("Reaction1",4,5,single).
bondMappingReactant("Reaction1",5,6,double).
bondMappingReactant("Reaction1",6,7,single).

atomMappingReactant("Reaction1",c,11).
atomMappingReactant("Reaction1",s,12).
atomMappingReactant("Reaction1",o,13).
atomMappingReactant("Reaction1",o,14).
bondMappingReactant("Reaction1",11,12,single).
bondMappingReactant("Reaction1",11,14,single).
bondMappingReactant("Reaction1",11,13,double).

atomMappingProduct("Reaction1",c,1).
atomMappingProduct("Reaction1",c,2).
atomMappingProduct("Reaction1",c,3).
atomMappingProduct("Reaction1",c,4).
atomMappingProduct("Reaction1",c,5).
atomMappingProduct("Reaction1",c,6).
atomMappingProduct("Reaction1",c,7).
atomMappingProduct("Reaction1",o,8).
atomMappingProduct("Reaction1",o,9).
atomMappingProduct("Reaction1",o,10).
atomMappingProduct("Reaction1",c,11).
atomMappingProduct("Reaction1",s,12).
atomMappingProduct("Reaction1",o,13).
atomMappingProduct("Reaction1",o,14).
bondMappingProduct("Reaction1",1,2,single).
bondMappingProduct("Reaction1",1,9,single).
bondMappingProduct("Reaction1",1,8,double).
bondMappingProduct("Reaction1",2,3,single).
bondMappingProduct("Reaction1",3,12,single).
bondMappingProduct("Reaction1",3,4,single).
bondMappingProduct("Reaction1",4,10,single).
bondMappingProduct("Reaction1",4,5,single).
bondMappingProduct("Reaction1",5,6,double).
bondMappingProduct("Reaction1",6,7,single).
bondMappingProduct("Reaction1",12,11,single).
bondMappingProduct("Reaction1",11,14,single).
bondMappingProduct("Reaction1",11,13,double).
