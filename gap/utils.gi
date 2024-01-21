#
# Author:   Suhas Vittal
# date:     19 January 2024
#
# Utility code for QEC.
#

LoadPackage("Gauss");

f_ListOrders :=
    function(G)
        local orders, x;
        orders := List(Elements(G), x -> Order(x));

        Print("Orders:");
        for x in orders do
            Print(" ", x);
        od;
        Print("\n");
    end;

# Returns true if H is suspected to be Torsion-Free.
f_IsTorsionFree := 
    function(H)
        local g, ord, result;
        for g in GeneratorsOfGroup(H) do
            result := IO_CallWithTimeoutList(rec(seconds := 20), Order, [g]);
            if result[1] then
                ord := result[2];
                if not (ord = infinity) then
                    return false;
                fi;
            fi;
        od;
        return true;
    end;

# Returns a list of vectors corresponding to each check using the coset
# representations.
f_ComputeChecks :=
    function(check_cosets, qubit_cosets)
        local F, G, H, common, v_check_list, v_check;

        F := GF(2);

        v_check_list := [];
        for G in check_cosets do
            v_check := [];
            for H in qubit_cosets do
                common := Intersection(AsList(H), AsList(G));
                if Length(common) > 0 then
                    Add(v_check, One(F));
                else
                    Add(v_check, Zero(F));
                fi;
            od;
            Add(v_check_list, v_check);
        od;
        return v_check_list;
    end;

f_Weight :=
    function(check)
        local n_nonzero, x;
        n_nonzero := 0;
        for x in check do
            if x = One(GF(2)) then
                n_nonzero := n_nonzero+1;
            fi;
        od;
        return n_nonzero;
    end;

# Returns a list of generators given the stabilizers.
f_ComputeStabilizerGenerators :=
    function(check_list)
        local M;
        # Perform a RREF to compute the stabilizers.
        M := EchelonMat(check_list);
        return check_list{ Filtered(M.heads, x -> x > 0) };
    end;

# Returns a record of X and Z operators (.x and .z respectively).
f_ComputeOperators :=
    function(x_check_list, z_check_list)
        local im_x,
                im_z,
                ker_x,
                ker_z,
                z_op_list,
                x_op_list;

        im_x := VectorSpace(GF(2), x_check_list);
        im_z := VectorSpace(GF(2), z_check_list);

        ker_x := VectorSpace(GF(2), NullspaceMat(TransposedMat(x_check_list)));
        ker_z := VectorSpace(GF(2), NullspaceMat(TransposedMat(z_check_list)));
        
        x_op_list := Filtered(BasisVectors(Basis(ker_x)), v -> not (v in im_x));
        z_op_list := Filtered(BasisVectors(Basis(ker_z)), v -> not (v in im_z));

        return rec(x := x_op_list, z := z_op_list);
    end;

f_AreCommuting :=
    function(v, w)
        return (v * w) = Zero(GF(2));
    end;

# Prints generators for tanner graph.
f_PrintGenerators :=
    function(stabilizers, type, tabs)
        local i, j, s, x;

        i := 0;
        for s in stabilizers do
            j := 0;
            Print(tabs, type, i);
            for x in s do
                if x = One(GF(2)) then
                    Print(",", j);
                fi;
                j := j+1;
            od;
            Print("\n");
            i := i+1;
        od;
    end;
