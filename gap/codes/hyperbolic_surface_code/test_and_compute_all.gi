#
# Author:   Suhas Vittal
# date:     19 January 2024 
#
# Code for generating Hyperbolic Surface Codes.
#
# Precondition: run init_surface.gi
#

verbose := true;

use_group := 1;
while use_group <= Length(normal_low_index_subgroups) do
    Print("Subgroup ", use_group, ":\n");

    Read("gap/codes/hyperbolic_surface_code/test_subgroup.gi");
    Read("gap/codes/hyperbolic_surface_code/compute_code.gi");

    use_group := use_group+1;
od;
