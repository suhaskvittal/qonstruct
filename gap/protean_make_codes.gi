# author: Suhas Vittal

# Surface Codes:
params := [
    [3, 8],
    [4, 5],
    [4, 6],
    [5, 5],
    [5, 6],
    [6, 6]
];

max_index := 4000;
for arr in params do
    r := arr[1];
    s := arr[2];
    Read("gap/codes/hyperbolic/surface_code.gi");
od;

# Color Codes:
params := [
    [3, 10],
    [3, 8],
    [4, 10],
    [4, 6],
    [4, 8],
    [5, 6]
];

max_index := 2000;
for arr in params do
    r := arr[1];
    s := arr[2];
    Read("gap/codes/hyperbolic/color_code.gi");
od;

