# Author: Suhas Vittal

i_ := E(4);
I := [[1,0,0,0],
      [0,1,0,0],
      [0,0,1,0],
      [0,0,0,1]];

KRON :=
    function(U)
        return [ KroneckerProduct(U, I), KroneckerProduct(I, U) ];        
    end;

# Single ququart gates:

H := ER(1/2) * [[1,1,0,0],
                [1,-1,0,0],
                [0,0,1,1],
                [0,0,1,-1]];
S := [[1,0,0,0],
      [0,E(4),0,0],
      [0,0,-E(4),0],
      [0,0,0,-1]];
L := [[0,1,0,0],
      [0,0,1,0],
      [0,0,0,1],
      [1,0,0,0]];
PX := [[0,1,0,0],
       [1,0,0,0],
       [0,0,0,E(4)],
       [0,0,E(4),0]];
PZ := [[1,0,0,0],
       [0,-1,0,0],
       [0,0,E(4),0],
       [0,0,0,-E(4)]];
PY := E(4) * PX * PZ;

CZ := 
       [
        [   1,0,0,0,    0,0,0,0,    0,0,0,0,    0,0,0,0     ],
        [   0,1,0,0,    0,0,0,0,    0,0,0,0,    0,0,0,0     ],
        [   0,0,1,0,    0,0,0,0,    0,0,0,0,    0,0,0,0     ],
        [   0,0,0,1,    0,0,0,0,    0,0,0,0,    0,0,0,0     ],

        [   0,0,0,0,    1,0,0,0,    0,0,0,0,    0,0,0,0     ],
        [   0,0,0,0,    0,-1,0,0,   0,0,0,0,    0,0,0,0     ],
        [   0,0,0,0,    0,0,1,0,    0,0,0,0,    0,0,0,0     ],
        [   0,0,0,0,    0,0,0,-1,   0,0,0,0,    0,0,0,0     ],

        [   0,0,0,0,    0,0,0,0,    1,0,0,0,    0,0,0,0     ],
        [   0,0,0,0,    0,0,0,0,    0,1,0,0,    0,0,0,0     ],
        [   0,0,0,0,    0,0,0,0,    0,0,1,0,    0,0,0,0     ],
        [   0,0,0,0,    0,0,0,0,    0,0,0,1,    0,0,0,0     ],

        [   0,0,0,0,    0,0,0,0,    0,0,0,0,    1,0,0,0     ],
        [   0,0,0,0,    0,0,0,0,    0,0,0,0,    0,-1,0,0    ],
        [   0,0,0,0,    0,0,0,0,    0,0,0,0,    0,0,1,0     ],
        [   0,0,0,0,    0,0,0,0,    0,0,0,0,    0,0,0,-1    ]
      ];

