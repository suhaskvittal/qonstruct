# Utility functions for files in this directory.

GetCycle :=
    function(x, offset)
        local arr, i;

        arr := [];
        i := 1;
        while i <= x do
            Add(arr, i+offset); 
            i := i+1;
        od;
        return CycleFromList(arr);
    end;

PermGroup := 
    function(arr)
        local _arr, i, off;

        _arr := [];
        i := 1;
        off := 0;
        while i <= Length(arr) do
            Add(_arr, GetCycle(arr[i], off));
            off := off+arr[i];
            i := i+1;
        od;
        return Group(_arr);
    end;

ComputeIndicatorVectors :=
    function(check_cosets, qubit_cosets)
        local G, H, common, v_check_list, v_check;

        v_check_list := [];
        for G in check_cosets do
            v_check := [];
            for H in qubit_cosets do
                common := Intersection(AsList(H), AsList(G));
                if Length(common) > 0 then
                    Add(v_check, One(GF(2))); else
                    Add(v_check, Zero(GF(2)));
                fi;
            od;
            Add(v_check_list, v_check);
        od;
        return v_check_list;
    end;

CodeFilename :=
    function(folder, n, k, dx, dz)
        local filename, d;
        if (dx < dz) then 
            d := dx; 
        else 
            d := dz;
        fi;
        filename := Concatenation(folder, "/", 
                                    String(n), "_", String(k), "_", String(d),
                                    ".txt");
        return filename;
    end;

GetFaceSet :=
    function(g, a, b, r, s)

    end;

FacesShareEdge :=
    function(f1, f2)
        return Length(Intersection(AsList(f1), AsList(f2))) > 0;
    end;

TilingIsBipartite :=
    function(faces)
        local colors, i, j, f1, f2, tries;

        # Colors tracks the color of the face. If -1, then the color is not set.
        colors := List(faces, x -> -1);
        colors[1] := 0;

        tries := 0;
        while ForAny(colors, x -> x = -1) and tries < 2 do
            i := 1;
            while i <= Length(faces) do
                # If the color is not set, then skip this color.
                if colors[i] = -1 then
                    i := i+1;
                    continue;
                fi;
                # We are essentially checking two-colorability.
                f1 := faces[i];
                j := i+1;
                while j <= Length(faces) do
                    f2 := faces[j];
                    if not FacesShareEdge(f1, f2) then
                        j := j+1;
                        continue;
                    fi;
                    if colors[j] = -1 then
                        colors[j] := 1 - colors[i];
                    fi;
                    if not (colors[i] = 1 - colors[j]) then
                        return false;
                    fi;
                    j := j+1;
                od;
                i := i+1;
            od;
            tries := tries+1;
        od;
        if tries = 2 then
            Print("\tWarning: could not converge to answer.\n");
            return false;
        fi;
        return true;
    end;
