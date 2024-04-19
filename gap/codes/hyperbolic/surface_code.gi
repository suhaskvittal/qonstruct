#
# Author:   Suhas Vittal
# date:     19 January 2024 
#
# Code for generating Hyperbolic Surface Codes.
#
# Load:
#   r and s -- the Schlaffi symbols for the surface.
#   max_index -- max index of normal subgroup

# Get functions.
Read("gap/utils/qec.gi");
Read("gap/utils/operations.gi");
Read("gap/codes/hyperbolic/base.gi");

LoadPackage("LINS");

# Compute normal subgroups.

#G_base := FreeGroup(3);
#G_base := G_base / [ G_base.1^2, G_base.2^2, G_base.3^2, 
#                   (G_base.1*G_base.2)^2, 
#                   (G_base.2*G_base.3)^r,
#                   (G_base.3*G_base.1)^s ];
#G_rs := Subgroup(G_base, [G_base.2*G_base.3, G_base.3*G_base.1]);
G := FreeGroup(2);
G := G / [ G.1^r, G.2^s, (G.1*G.2)^2 ];

gr := LowIndexNormalSubgroupsSearch(G, max_index);
lins := ComputedNormalSubgroups(gr);

Print("# of normal subgroups = ", Length(lins), "\n");

G := Grp(LinsRoot(gr));

# Make folder for codes.
code_folder := Concatenation("codes/hysc/", String(r), "_", String(s));
IO_mkdir(code_folder, 448);

# Iterate through subgroups.
i := 2;
visited := [];
while i <= Length(lins) do
    if Index(lins[i]) in visited then
        i := i+1;
        continue;
    fi;
    Print("Index: ", Index(lins[i]), "\n");
    AddSet(visited, Index(lins[i]));

    H := Grp(lins[i]);
    # Compute quotient group.
    G_mod_H := G / H;
    gens := GeneratorsOfGroup(G_mod_H);
    # Print generator orders.
    rho := gens[1];
    sig := gens[2];
    tau := (rho*sig)^-1;
    Print("Orders of tau, sigma, and rho: ", Order(tau), ", ", Order(sig), ", ", Order(rho), "\n");
    if not (Order(tau) = 2 and Order(sig) = s and Order(rho) = r) then
        i := i+1;
        continue;
    fi;
    # Get isomorphic subgroups:
    H_bc := Subgroup(G_mod_H, [rho]);
    H_ca := Subgroup(G_mod_H, [sig]);
    H_ab := Subgroup(G_mod_H, [tau]);
    # Read code construction file.
    faces := RightCosets(G_mod_H, H_bc);
    vertices := RightCosets(G_mod_H, H_ca);
    edges := RightCosets(G_mod_H, H_ab);

    n_data := Length(edges);
    genus := 2 - (Length(vertices) - Length(edges) + Length(faces));

    # Compute checks. A data qubit is in the support of an X (Z) check if
    # its edge coset shares an element with the check's face (vertex) coset.
    x_checks := ComputeIndicatorVectorsSC(vertices, edges);
    z_checks := ComputeIndicatorVectorsSC(faces, edges);

    # Compute stabilizer generators:
    x_checks := ComputeStabilizerGenerators(x_checks);
    z_checks := ComputeStabilizerGenerators(z_checks);

    # Compute logical operators:
    x_operators := CssLXZOperators(x_checks, z_checks);
    z_operators := CssLXZOperators(z_checks, x_checks);

    n_x_checks := Length(x_checks);
    n_z_checks := Length(z_checks);
    n_checks := n_x_checks + n_z_checks;
    n_x_ops := Length(x_operators);
    n_z_ops := Length(z_operators);

    dx := MinOperatorWeight(x_operators);
    dz := MinOperatorWeight(z_operators);

    # Check if the code is valid. If so, write it to a file.
    if n_data - n_checks = genus
    and n_x_ops = genus
    and n_z_ops = genus
    then
        Print("\tData qubits: ", n_data, "\n");
        Print("\tLogical qubits: ", genus, "\n");
        Print("\tChecks: ", n_checks, "\n");
        Print("\t\tX: ", n_x_checks, ", Z: ", n_z_checks, "\n");
        Print("\tX/Z Operators: ", Length(x_operators), ", ", Length(z_operators), "\n");
        Print("\t\tX min weight: ", dx, "\n");
        Print("\t\tZ min weight: ", dz, "\n");

        WriteCodeToFile(CodeFilename(code_folder, n_data, genus),    
                            x_checks,
                            z_checks,
                            x_operators,
                            z_operators);
        done := true;
    fi;
    i := i+1;
od;
