#
# Author:   Suhas Vittal
# date:     19 January 2024 
#
# Code for generating Hyperbolic Surface Codes.
#
# Precondition: run test_subgroup.gi.
#
# Load: verbose (true if print out debugging messages).

f_ComputeIndicatorVectors :=
    function(check_cosets, qubit_cosets)
        local G, H, common, v_check_list, v_check;

        v_check_list := [];
        for G in check_cosets do
            v_check := [];
            for H in qubit_cosets do
                common := Intersection(AsList(H), AsList(G));
                if Length(common) > 0 then
                    Add(v_check, One(GF(2)));
                else
                    Add(v_check, Zero(GF(2)));
                fi;
            od;
            Add(v_check_list, v_check);
        od;
        return v_check_list;
    end;

done := false;

ii := 1;
while ii <= Length(iso_bc) do
    jj := 1;
    while jj <= Length(iso_ca) do
        kk := 1;
        while kk <= Length(iso_ab) do
            H_bc := Image(iso_bc[ii]);
            H_ca := Image(iso_ca[jj]);
            H_ab := Image(iso_ab[kk]);

            faces := RightCosets(G_rs_mod_H, H_bc);
            vertices := RightCosets(G_rs_mod_H, H_ca);
            edges := RightCosets(G_rs_mod_H, H_ab);

            genus := 2 - (Length(vertices) - Length(edges) + Length(faces));

            n_data := Length(edges);
            # Compute checks. A data qubit is in the support of an X (Z) check if
            # its edge coset shares an element with the check's face (vertex) coset.
            x_checks := f_ComputeIndicatorVectors(vertices, edges);
            z_checks := f_ComputeIndicatorVectors(faces, edges);
            # Compute stabilizer generators:
            x_checks := f_ComputeStabilizerGenerators(x_checks);
            z_checks := f_ComputeStabilizerGenerators(z_checks);

            n_x_checks := Length(x_checks);
            n_z_checks := Length(z_checks);
            n_checks := n_x_checks + n_z_checks;

            invalid := not (n_data - n_checks = genus);
            if invalid then
                if verbose then
                    Print("\tFailed as ", genus, " != ", n_data, " - ", n_checks, 
                            "(", n_data-n_checks, ")\n");
                fi;
                kk := kk+1;
            else 
                Print("\tData qubits: ", n_data, "\n");
                Print("\tLogical qubits: ", genus, "\n");
                Print("\tChecks: ", n_checks, "\n");
                Print("\t\tX: ", n_x_checks, ", Z: ", n_z_checks, "\n");
                # Now, split the checks into X and Z checks. We will now try to
                # compute the logical operators.
                x_operators := f_CssLXZOperators(x_checks, z_checks);
                z_operators := f_CssLXZOperators(z_checks, x_checks);

                break;
            fi;
        od;
        if done then break; fi;
        jj := jj+1;
    od;
    if done then break; fi;
    ii := ii+1;
od;
