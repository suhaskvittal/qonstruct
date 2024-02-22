#
# Author:   Suhas Vittal
# date:     19 January 2024
#
# Utility code for QEC.
#

LoadPackage("IO");

OperatorWeight :=
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

ComputeStabilizerGenerators :=
    function(check_list)
        local M;
        M := SemiEchelonMat(check_list);
        return check_list{ Filtered(M.heads, x -> x > 0) };
    end;

AreCommuting :=
    function(v, w)
        return (v * w) = Zero(GF(2));
    end;

# Computes the logical X/Z operators. To compute logical X operators,
# the first argument should be the X checks and the second should be
# the Z checks. For logical Z, vice versa.
CssLXZOperators := 
    function(same_type_checks, opposing_checks)
        local im_osame, im_oopp, ker, op_list;

        im_osame := VectorSpace(GF(2), same_type_checks); 
        
        op_list := [];
        Append(op_list, same_type_checks);
        Append(op_list, NullspaceMat(TransposedMat(opposing_checks)));

        ker := VectorSpace(GF(2), op_list);
        return Filtered(BasisVectors(Basis(ker)), v -> not (v in im_osame));
    end;

MinOperatorWeight :=
    function(operators) 
        local min_w, op, w;

        min_w := OperatorWeight(operators[1]);
        for op in operators do
            w := OperatorWeight(op);
            if w < min_w then   min_w := w; fi;
        od;
        return min_w;
    end;

# Prints generators for tanner graph.
PrintGenerators :=
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

WriteChecksAndOperatorsToFile :=
    function(fout, checks, operators, type)
        local line, i, j, s, x;

        i := 0;
        for s in checks do
            line := Concatenation(type, String(i));
            j := 0;
            for x in s do
                if x = One(GF(2)) then
                    line := Concatenation(line, ",", String(j));
                fi;
                j := j+1;
            od;
            i := i+1;

            IO_WriteLine(fout, line);
        od;
        i := 0;
        for s in operators do
            line := Concatenation("O", type, String(i));
            j := 0;
            for x in s do
                if x = One(GF(2)) then
                    line := Concatenation(line, ",", String(j));
                fi;
                j := j+1;
            od;
            i := i+1;

            IO_WriteLine(fout, line);
        od;
    end;

WriteCodeToFile :=
    function(file, x_checks, z_checks, x_operators, z_operators)
        local fout;

        fout := IO_File(file, "w");
        Print("Creating file ", file, ": ", IsFile(fout), "\n");
        WriteChecksAndOperatorsToFile(fout, x_checks, x_operators, "X");
        WriteChecksAndOperatorsToFile(fout, z_checks, z_operators, "Z");
    end;
