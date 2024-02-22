# Author: Suhas Vittal

Read("gap/gates.gi");

gens := [CX];
Append(gens, KRON(H));
Append(gens, KRON(S));
Append(gens, KRON(L));
 
# Search for normal subgroup.
#CG := GroupByGenerators(gens);
CG := Group(H, S, L);
nsgs := NormalSubgroups(CG);
# Find smallest normal subgroup with X, Y, and Z.
PX := [[0,1,0,0],
       [1,0,0,0],
       [0,0,0,1],
       [0,0,1,0]];
PZ := [[1,0,0,0],
       [0,-1,0,0],
       [0,0,1,0],
       [0,0,0,-1]];

k := 1;
min_size := -1;
while k <= Length(nsgs) do
    SG := nsgs[k];
    n := Length(Elements(SG));
    ng := Length(GeneratorsOfGroup(SG));
    if min_size > 0 and n > min_size then
        k := k+1;
        continue;
    fi;

    is_good := (S in SG) and (L in SG) and (PX in SG) and (PZ in SG);
    if is_good and (min_size < 0 or n < min_size) then
        Print("Found group with ", n, " elements and ", ng, " generators.\n");
        min_size := n;
        best_grp := SG;
    fi;
    k := k+1;
od;
