# Code for generating Hyperbolic Surface Codes.
#
# Dependency: LINS (https://github.com/gap-packages/LINS)
# 
# Precondition: load r and s beforehand.

LoadPackage("LINS");

f_IsTorsionFree := 
    function(H)
        local generators, g, ord, result;
        Print("\tchecking if torsion free...\n");

        generators := GeneratorsOfGroup(H);
        for g in generators do
            result := IO_CallWithTimeoutList(rec(seconds := 10), Order, [g]);
            if result[1] then
                ord := result[2];
                if not (ord = infinity) then
                    Print("\tfound generator of order ", ord, "\n");
                    return false;
                fi;
            fi;
        od;
        return true;
    end;

f_LeftCoset := 
    function(x, H)
        local coset_inv, coset, y;
        # Get right coset of inverse (that is, Hx^-1). Convert this to left
        # coset.
        coset_inv := RightCoset(H, x^-1);
        coset := [];
        for y in coset_inv do
            AddSet(coset, y^-1);
        od;
        return coset;
    end;
    

r := 4;
s := 5;

G_base := FreeGroup(3);
G_base := G_base / [ G_base.1^2, G_base.2^2, G_base.3^2, 
                    (G_base.1*G_base.2)^2, 
                    (G_base.2*G_base.3)^r,
                    (G_base.3*G_base.1)^s ];
gA := G_base.1;
gB := G_base.2;
gC := G_base.3;

gAB := gA*gB;
gBC := gB*gC;
gCA := gC*gA;

G_rs := Subgroup(G_base, [G_base.2*G_base.3, G_base.3*G_base.1]);

# Search for torsion-free normal subgroups.

tfn_subgroup_found := false;
tfn_subgroup := 0;

gr := LowIndexNormalSubgroupsSearch(G_rs, 40);
normal_low_index_subgroups := ComputedNormalSubgroups(gr);

# A quirk of LINS: the root of the graph is different from G_rs,
# so we need to update G_rs.

gr_root := LinsRoot(gr);
_G_rs := G_rs;
G_rs := Grp(gr_root);

i := 0;
for H in normal_low_index_subgroups do
    H_index := Index(H);
    Print("Subgroup ", i, ", index: ", H_index, "\n");
    i := i+1;

    _H := Grp(H);
    # Need to convert to subgroup of G_rs.
    H := SubgroupNC(G_rs, GeneratorsOfGroup(_H));
    # Check if H is Torsion Free.
    if f_IsTorsionFree(H) then
        tfn_subgroup_found := true;
        tfn_subgroup := H;
    fi;
od;

if tfn_subgroup_found then
    G_rs_mod_H := G_rs / tfn_subgroup;

    H_bc := Subgroup(G_rs_mod_H, [G_rs.1]);
    H_ca := Subgroup(G_rs_mod_H, [G_rs.2]);
    H_ab := Subgroup(G_rs_mod_H, [(G_rs.1*G_rs.2)^-1]);

    faces := [];
    vertices := [];
    edges := [];

    for x in G_rs_mod_H do
        for y in f_LeftCoset(x, H_bc) do    AddSet(faces, y);       od;
        for y in f_LeftCoset(x, H_ca) do    AddSet(vertices, y);    od;
        for y in f_LeftCoset(x, H_ab) do    AddSet(edges, y);       od;
    od;

    genus := 2 - (Length(vertices) - Length(edges) + Length(faces));
    Print("# vertices, edges, faces = ", Length(vertices), Length(edges), Length(faces), "\n");
    Print("Logical qubits: ", genus, "\n");
fi;
