# Utility functions for files in this directory.

ComputeIndicatorVectorsSC :=
    function(check_cosets, qubit_cosets)
        local G, H, v_check_list, v_check;

        v_check_list := [];
        for G in check_cosets do
            v_check := [];
            for H in qubit_cosets do
                if Length(Intersection(AsList(H), AsList(G))) > 0 then
                    Add(v_check, One(GF(2)));
                else
                    Add(v_check, Zero(GF(2)));
                fi;
            od;
            Add(v_check_list, v_check);
        od;
        return v_check_list;
    end;

ComputeIndicatorVectorsCC :=
    function(G, faces)
        local check_list, check, f, x;

        check_list := [];
        for f in faces do 
            check := [];
            for x in Elements(G) do
                if x in f then
                    Add(check, One(GF(2)));
                else
                    Add(check, Zero(GF(2)));
                fi;
            od;
            Add(check_list, check);
        od;
        return check_list;
    end;

CodeFilename :=
    function(folder, n, k)
        local filename;
        filename := Concatenation(folder, "/", 
                                    String(n), "_", String(k), ".txt");
        return filename;
    end;

Rot1 :=
    function(b)
        local arr, i;
        arr := [];
        i := 0;
        while i < Order(b) do
            AddSet(arr, b^i);
            i := i+1;
        od;
        return arr;
    end;

Rot2 :=
    function(a, b)
        local arr, i, k;
        arr := [];
        i := 0;
        k := 0;
        while i < 2*Order(a*b) do
            if i mod 2 = 0 then
                AddSet(arr, (a*b)^k);
                k := k+1;
            else
                AddSet(arr, (a*b)^k * a);
            fi;
            i := i+1;
        od;
        return arr;
    end;

GetCosets :=
    function(G, rots)
        local arr, _arr, x, y;
        arr := [];
        for x in Elements(G) do
            _arr := [];
            for y in rots do
                AddSet(_arr, x*y);
            od;
            AddSet(arr, _arr);
        od;
        return arr;
    end;

FacesShareEdge :=
    function(f1, f2)
        return Length(Intersection(f1, f2)) = 2;
    end;

TilingIsBipartite :=
    function(faces)
        local colors,
                dfs,
                visited,
                i,
                j,
                f1,
                f2;

        # Colors tracks the color of the face. If -1, then the color is not set.
        colors := List(faces, x -> -1);
        colors[1] := 0;

        dfs := [1];
        visited := [];
        while Size(dfs) > 0 do
            i := dfs[Size(dfs)];
            Remove(dfs);
            if i in visited then
                continue;
            fi;
            f1 := faces[i];
            # Iterate through all the faces. Only visit those adjacent to f1.
            j := 1;
            while j <= Size(faces) do
                if i = j then 
                    j := j+1;
                    continue;
                fi;
                f2 := faces[j];
                if not FacesShareEdge(f1, f2) then
                    j := j+1;
                    continue;
                fi;
                if colors[j] >= 0 then
                    if colors[i] = colors[j] then
                        return false;
                    else
                        j := j+1;
                        continue;
                    fi;
                fi;
                colors[j] := 1-colors[i];
                Add(dfs, j);
            od;
            AddSet(visited, i);
        od;
        return true;
    end;

