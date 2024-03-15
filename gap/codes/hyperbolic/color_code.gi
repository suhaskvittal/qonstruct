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

# Get functions.
Read("gap/utils/qec.gi");
Read("gap/utils/operations.gi");
Read("gap/codes/hyperbolic/base.gi");

LoadPackage("LINS");

# Construct tiling:
G := FreeGroup(2);
G := G / [ G.1^r, G.2^s, (G.1*G.2)^2 ];

# Get low-index normal subgroup of G. We will use this to make a quotient group.
gr := LowIndexNormalSubgroupsSearch(G, max_index);
lins := ComputedNormalSubgroups(gr);

Print("# of normal subgroups = ", Length(lins), "\n");

G := Grp(LinsRoot(gr));

code_folder := Concatenation("codes/hycc/", String(r), "_", String(s));
IO_mkdir(code_folder, 448);

# Iterate through subgroups.
i := 1;
visited := [];
while i <= Length(lins) do
    n_data := Index(lins[i]); # This is the number of qubits.

    if Index(lins[i]) in visited then
        i := i+1;
        continue;
    fi;
    Print("Index: ", Index(lins[i]), "\n");
    AddSet(visited, Index(lins[i]));

    H := Grp(lins[i]);

    G_mod_H := G / H;
    gens := GeneratorsOfGroup(G_mod_H);
    # Generators of the finite group:
    a := gens[1];
    b := gens[2];
    ab := a*b;
    
    Print("\tGenerator orders:", Order(a), ", ", Order(b), ", ", Order(ab), "\n");

    if not (Order(a) = r and Order(b) = s and Order(ab) = 2) then
        i := i+1;
        continue;
    fi;

    _f1 := Rot1(b);
    _f2 := Rot2(ab, b);
    _f3 := Rot2(ab, b^-1);

    f1 := GetCosets(G_mod_H, _f1);
    f2 := GetCosets(G_mod_H, _f2);
    f3 := GetCosets(G_mod_H, _f3);

    f23 := f2;
    Append(f23, f3);
    is_bp := TilingIsBipartite(f23);
    if not is_bp then
        i := i+1;
        continue;
    fi;
    # Identify each face (in f1, f2, f3) as a check operator, and each element
    # of the coset is a qubit.
    all_faces := f1;
    Append(all_faces, f23);

    plaquettes := ComputeIndicatorVectorsCC(G_mod_H, all_faces);
    plaquettes := ComputeStabilizerGenerators(plaquettes);

    if n_data <= 2*Length(plaquettes) then
        i := i+1;
        continue;
    fi;

    Print("\tPlaquettes = ", Length(plaquettes), "\n");

    # Compute logical operators of the code.
    operators := CssLXZOperators(plaquettes, plaquettes);
    
    n_checks := 2*Length(plaquettes);
    n_ops := Length(operators);
    n_log := n_data - n_checks;
    d := MinOperatorWeight(operators);

    Print("\tData qubits: ", n_data, "\n");
    Print("\tLogical qubits: ", n_log, "\n");
    Print("\tChecks: ", n_checks, "\n");
    Print("\tOperators: ", n_ops, "\n");
    Print("\t\tMin weight: ", d, "\n");

    WriteCodeToFile(CodeFilename(code_folder, n_data, n_log),    
                        plaquettes,
                        plaquettes,
                        operators,
                        operators);

    i := i+1;
od;
