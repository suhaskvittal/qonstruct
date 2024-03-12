# Author: Suhas Vittal

Read("gap/leakage_sim/gates.gi");

# Search for normal subgroup.
CG := Group(H, S, L, -I);

#CG2 := GroupByGenerators(gens);

# It's very hard to get normal subgroups of CG2.
nsgs := NormalSubgroups(CG);
# Find smallest normal subgroup with X, Y, and Z.

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

    is_good := (H*S*H in SG) and (S in SG);
    if is_good and (min_size < 0 or n < min_size) then
        Print("Found group with ", n, " elements and ", ng, " generators.\n");
        min_size := n;
        best_grp := SG;
    fi;
    k := k+1;
od;
