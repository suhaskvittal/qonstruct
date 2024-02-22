#
# Author:   Suhas Vittal
# date:     19 January 2024 
#
# Code for generating Hyperbolic Color Codes.
#
# Load:
#   r and s -- the Schlaffi symbols for the surface.
#   max_index -- max index of normal subgroup
#

r := 3;
s := 8;
max_index := 2000;

# Get functions.
Read("gap/utils/qec.gi");
Read("gap/utils/operations.gi");
Read("gap/codes/hyperbolic/base.gi");

LoadPackage("LINS");

# Construct tiling:
G := FreeGroup(2);
G := G / [ G.1^2, G.2^r, (G.1*G.2)^s ];
#G_base := FreeGroup(3);
#G_base := G_base / [ G_base.1^2, G_base.2^2, G_base.3^2, 
#                   (G_base.1*G_base.2)^2, 
#                   (G_base.2*G_base.3)^r,
#                   (G_base.3*G_base.1)^s ];
#G := Subgroup(G_base, [G_base.2*G_base.3, G_base.3*G_base.1]);

# Get low-index normal subgroup of G. We will use this to make a quotient group.
gr := LowIndexNormalSubgroupsSearch(G, max_index);
lins := ComputedNormalSubgroups(gr);

Print("# of normal subgroups = ", Length(lins), "\n");

G := Grp(LinsRoot(gr));

# Iterate through subgroups.
i := 1;
while i <= Length(lins) do
    Print("Index: ", Index(lins[i]), "\n");
    H := Grp(lins[i]);

    G_mod_H := G / H;
    # Get isomorphic subgroups of interest:
    iso_1 := IsomorphicSubgroups(G_mod_H, PermGroup([r]));
    iso_2 := IsomorphicSubgroups(G_mod_H, PermGroup([2, r]));
    iso_3 := IsomorphicSubgroups(G_mod_H, PermGroup([2, s]));
    Print("Subgroups: ", Length(iso_1), ", ", Length(iso_2), ", ", Length(iso_3), "\n");
    # Iterate through all isomorphisms.
    done := false;
    ii := 1;
    while ii <= Length(iso_1) do
        jj := 1;
        while jj <= Length(iso_2) do
            kk := 1;
            while kk <= Length(iso_3) do
                # Get subgroups.
                H1 := Image(iso_1[ii]);
                H2 := Image(iso_2[jj]);
                H3 := Image(iso_3[kk]);
                faces1 := RightCosets(G_mod_H, H1);
                faces2 := RightCosets(G_mod_H, H2);
                faces3 := RightCosets(G_mod_H, H3);

                faces23 := Union(faces2, faces3);
                Print(Length(faces23), ", ", Length(faces2), ", ", Length(faces3), "\n");
                if TilingIsBipartite(faces23) then
                    Print("\tTiling is Bipartite.\n");
                    done := true;
                else
                    Print("\tTiling is not bipartite.\n");
                fi;
                
                if done then break; fi;
                kk := kk+1;
            od;
            if done then break; fi;
            jj := jj+1;
        od;
        if done then break; fi;
        ii := ii+1;
    od;
    i := i+1;
od;
