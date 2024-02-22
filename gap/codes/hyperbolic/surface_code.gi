#
# Author:   Suhas Vittal
# date:     19 January 2024 
#
# Code for generating Hyperbolic Surface Codes.
#
# Load:
#   r and s -- the Schlaffi symbols for the surface.
#   max_index -- max index of normal subgroup

r := 4;
s := 5;
max_index := 1000;

# Get functions.
Read("gap/utils/qec.gi");
Read("gap/utils/operations.gi");
Read("gap/codes/hyperbolic/base.gi");

LoadPackage("LINS");

# Compute normal subgroups.

G_base := FreeGroup(3);
G_base := G_base / [ G_base.1^2, G_base.2^2, G_base.3^2, 
                    (G_base.1*G_base.2)^2, 
                    (G_base.2*G_base.3)^r,
                    (G_base.3*G_base.1)^s ];
G_rs := Subgroup(G_base, [G_base.2*G_base.3, G_base.3*G_base.1]);

gr := LowIndexNormalSubgroupsSearch(G_rs, max_index);
lins := ComputedNormalSubgroups(gr);

Print("# of normal subgroups = ", Length(lins), "\n");

G_rs := Grp(LinsRoot(gr));

# Make folder for codes.
code_folder := Concatenation("codes/hysc/", String(r), "_", String(s));
IO_mkdir(code_folder, 448);

# Iterate through subgroups.
i := 1;
while i <= Length(lins) do
    Print("Index: ", Index(lins[i]), "\n");
    H := Grp(lins[i]);
    # Compute quotient group.
    G_rs_mod_H := G_rs / H;
    # Get isomorphic subgroups:
    iso_1 := IsomorphicSubgroups(G_rs_mod_H, PermGroup([r]));
    iso_2 := IsomorphicSubgroups(G_rs_mod_H, PermGroup([s]));
    iso_ab := IsomorphicSubgroups(G_rs_mod_H, Group((1,2)));
    Print("Subgroups: ", Length(iso_1), ", ", Length(iso_2), ", ", Length(iso_ab), "\n");
    # Iterate through all isomorphisms.
    done := false;
    ii := 1;
    while ii <= Length(iso_1) do
        jj := 1;
        while jj <= Length(iso_2) do
            kk := 1;
            while kk <= Length(iso_ab) do
                # Get subgroups.
                H_bc := Image(iso_1[ii]);
                H_ca := Image(iso_2[jj]);
                H_ab := Image(iso_ab[kk]);
                # Read code construction file.
                faces := RightCosets(G_rs_mod_H, H_bc);
                vertices := RightCosets(G_rs_mod_H, H_ca);
                edges := RightCosets(G_rs_mod_H, H_ab);

                n_data := Length(edges);
                genus := 2 - (Length(vertices) - Length(edges) + Length(faces));

                # Compute checks. A data qubit is in the support of an X (Z) check if
                # its edge coset shares an element with the check's face (vertex) coset.
                x_checks := ComputeIndicatorVectors(vertices, edges);
                z_checks := ComputeIndicatorVectors(faces, edges);

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
                and not (dx <= 2 or dz <= 2)
                then
                    Print("\tData qubits: ", n_data, "\n");
                    Print("\tLogical qubits: ", genus, "\n");
                    Print("\tChecks: ", n_checks, "\n");
                    Print("\t\tX: ", n_x_checks, ", Z: ", n_z_checks, "\n");
                    Print("\tX/Z Operators: ", Length(x_operators), ", ", Length(z_operators), "\n");
                    Print("\t\tX min weight: ", dx, "\n");
                    Print("\t\tZ min weight: ", dz, "\n");

                    WriteCodeToFile(CodeFilename(code_folder, n_data, genus, dx, dz),    
                                        x_checks,
                                        z_checks,
                                        x_operators,
                                        z_operators);
                    done := true;
                fi;
                if (done) then break; fi;
                kk := kk+1;
            od;
            if (done) then break; fi;
            jj := jj+1;
        od;
        if (done) then break; fi;
        ii := ii+1;
    od;
    i := i+1;
od;
