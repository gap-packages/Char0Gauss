InstallGlobalFunction(C0GAUSS_FoldList2,
    function(list, func, op, def)
    local k, s, old_s, r, i, len, n, nh, res, r1, r2;


    len := Length(list);
    # FIXME: We don't know a default value to return here.
    if len = 0 then
        return def;
    elif len = 1 then
        return op(def, func(list[1]));
    fi;

    res := List(list, func);
    k := len;
    s := 1;
    while k > 1 do
        r := k mod 2;
        old_s := s;
        k := QuoInt(k, 2);
        s := s * 2;
        i := s;
        while i <= k * s do
            if IsBound(res[i-old_s]) then
                r1 := res[i-old_s];
            else
                r1 := 1;
            fi;
            if IsBound(res[i]) then
                r2 := res[i];
            else
                r2 := 1;
            fi;
            res[i] := op(r1, r2);
            res[i-old_s] := 0;
            i := i + s;
        od;
        if r = 1 then
            k := k + 1;
            res[i] := res[i-old_s];
        fi;
    od;
    return op(def, res[ k * s ]);
end );

InstallGlobalFunction(MatIntFFESymm,
                      mat -> List(mat, IntFFESymm) );

InstallOtherMethod(IntFFESymm,
                   "for matrices",
                   [IsFFECollColl],
                   MatIntFFESymm );

InstallMethod( C0GAUSS_LeastCommonDenominator, "for a sparse matrix",
               [ IsSparseMatrixRep ],
function(mat)
    if RingOfDefinition(mat) <> Rationals then
        TryNextMethod();
    fi;
    return C0GAUSS_FoldList2(List(mat!.entries, r -> C0GAUSS_FoldList2(r, DenominatorRat, LcmInt, 1)),
				             IdFunc, LcmInt, 1);
end);

InstallMethod( C0GAUSS_LeastCommonDenominator, "for a sparse matrix",
               [ IsSparseMatrixRep, IsPosInt ],
function(mat, init)
    if RingOfDefinition(mat) <> Rationals then
        TryNextMethod();
    fi;
    return C0GAUSS_FoldList2(List(mat!.entries, r -> C0GAUSS_FoldList2(r, DenominatorRat, LcmInt, init)),
                            IdFunc, LcmInt, init);
end);

InstallMethod( C0GAUSS_LeastCommonDenominator, "for a sparse matrix",
               [ IsMatrix ],
function(mat)
    if DefaultFieldOfMatrix(mat) <> Rationals then
        TryNextMethod();
    fi;
    return C0GAUSS_FoldList2(List(mat, r -> C0GAUSS_FoldList2(r, DenominatorRat, LcmInt, 1)),
                            IdFunc, LcmInt, 1);
end);

InstallMethod( C0GAUSS_LeastCommonDenominator, "for a list-of-lists matrix",
                [ IsMatrix, IsPosInt ],
function(mat, init)
    if DefaultFieldOfMatrix(mat) <> Rationals then
        TryNextMethod();
    fi;
    return C0GAUSS_FoldList2(List(mat, r -> C0GAUSS_FoldList2(r, DenominatorRat, LcmInt, init)),
          IdFunc, LcmInt, init);
end);


InstallMethod( C0GAUSS_LeastCommonDenominator, "for a sparse matrix",
               [ IsSparseMatrixRep, IsPosInt ],
function(mat, init)
    if RingOfDefinition(mat) <> Rationals then
        TryNextMethod();
    fi;
    return C0GAUSS_FoldList2(List(mat!.entries, r -> C0GAUSS_FoldList2(r, DenominatorRat, LcmInt, init)),
                            IdFunc, LcmInt, init);
end);
