#
# Author:   Suhas Vittal
# date:     19 January 2024 
#
# Code for generating Hyperbolic Surface Codes.
#
# Precondition: run test_subgroup.gi.
#
# Load: verbose (true if print out debugging messages).

done := false;

ii := 1;
while ii <= Length(iso_bc) do
    jj := 1;
    while jj <= Length(iso_ca) do
        kk := 1;
        while kk <= Length(iso_ab) do
            if verbose then
                Print("Trying i, j, k = ", ii, ", ", jj, ", ", kk, "\n");
            fi;
            H_bc := Image(iso_bc[ii]);
            H_ca := Image(iso_ca[jj]);
            H_ab := Image(iso_ab[kk]);

            faces := RightCosets(G_rs_mod_H, H_bc);
            vertices := RightCosets(G_rs_mod_H, H_ca);
            edges := RightCosets(G_rs_mod_H, H_ab);

            genus := 2 - (Length(vertices) - Length(edges) + Length(faces));

            # Compute checks. A data qubit is in the support of an X (Z) check if
            # its edge coset shares an element with the check's face (vertex) coset.
            x_checks := f_ComputeChecks(vertices, edges);
            z_checks := f_ComputeChecks(faces, edges);
            # Validate checks.
            invalid := false;
            for check in x_checks do
                invalid := invalid or not (f_Weight(check) = s);
            od;
            for check in z_checks do
                invalid := invalid or not (f_Weight(check) = r);
            od;
            # Compute generators.
            x_checks := f_ComputeStabilizerGenerators(x_checks);
            z_checks := f_ComputeStabilizerGenerators(z_checks);

            n_data := Length(edges);
            n_x_checks := Length(x_checks);
            n_z_checks := Length(z_checks);

            invalid := invalid or not (n_data - n_x_checks - n_z_checks = genus);

            if invalid then
                if verbose then
                    Print("\tFailed as ", genus, " != ", n_data, " - ", n_x_checks, " - ", n_z_checks, 
                            "(", n_data-n_x_checks-n_z_checks, ")\n");
                fi;
                kk := kk+1;
            else 
                Print("\tData qubits: ", n_data, "\n");
                Print("\tLogical qubits: ", genus, "\n");
                Print("\tX, Z checks = ", n_x_checks, ", ", n_z_checks, "\n");
#               f_PrintGenerators(x_checks, "x", "\t\t");
#               f_PrintGenerators(z_checks, "z", "\t\t");
                vs := f_ComputeOperators(x_checks, z_checks);
                n_x_ops := Length(vs.x);
                n_z_ops := Length(vs.z);
                
                Print("\t# of X, Z operators = ", n_x_ops, ", ", n_z_ops, "\n");
                kk := kk+1;
                break;
            fi;
        od;
        if done then break; fi;
        jj := jj+1;
    od;
    if done then break; fi;
    ii := ii+1;
od;
