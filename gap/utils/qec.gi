#
# Author:   Suhas Vittal
# date:     19 January 2024
#
# Utility code for QEC.
#

f_OperatorWeight :=
    function(operator)
        local w, x;
        w := 0;
        for x in operator do
            if x = One(GF(2)) then
                w := w+1;
            fi;
        od;
        return w;
    end;

f_ComputeStabilizerGenerators :=
    function(check_list)
        local M;
        M := SemiEchelonMat(check_list);
        return check_list{ Filtered(M.heads, x -> x > 0) };
    end;

f_AreCommuting :=
    function(v, w)
        return (v * w) = Zero(GF(2));
    end;

# Computes the logical X/Z operators. To compute logical X operators,
# the first argument should be the X checks and the second should be
# the Z checks. For logical Z, vice versa.
f_CssLXZOperators := 
    function(same_type_checks, opposing_checks)
        local im_of_same, ker_of_opp, op_list;

        im_of_same := VectorSpace(GF(2), same_type_checks); 
        ker_of_opp := VectorSpace(GF(2), NullspaceMat(TransposedMat(opposing_checks)));

        op_list := Filtered(BasisVectors(Basis(ker_of_opp)), v -> not (v in im_of_same));
        op_list := Filtered(op_list, x -> f_OperatorWeight(x) > 0);
        return op_list;
    end;

# Prints generators for tanner graph.
f_PrintGenerators :=
    function(stabilizers, type, tabs)
        local i, j, s, x;

        i := 0;
        for s in stabilizers do
            Print(tabs, type, i);
            j := 0;
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
