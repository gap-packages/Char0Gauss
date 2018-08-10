#
# Solving linear equations over the integers/rationals by Dixon/Hensel lifting
#
# This is a GAP prototype which already works quite a bit faster than any code
# inside GAP, for reasonably large systems (thousands of equations)
#
InstallGlobalFunction(C0GAUSS_SetupMatVecsSystem_Padic,
function(mat, vecs, p, precision)
    local system, mmults, vmults, lcm, t, r, n;
    system := rec( mat := mat
                 , vecs := vecs
                 , number_variables := Length(mat[1])
                 , number_equations := Length(mat) );

    #  MakeIntSystem(system);
    Info(InfoChar0GaussLinearEq, 5,
         "MakeIntSystem2: computing denominator lcms" );

    t := NanosecondsSinceEpoch();
    mmults := List(system.mat, x -> C0GAUSS_FoldList2(x, DenominatorRat, LcmInt));
    vmults := List(system.vecs, x -> C0GAUSS_FoldList2(x, DenominatorRat, LcmInt));
    lcm := C0GAUSS_FoldList2(Concatenation(mmults, vmults), IdFunc, LcmInt);
    t := NanosecondsSinceEpoch() - t;
    Info(InfoChar0GaussLinearEq, 1, "computing the LCM took: ", t / 1000000., " msec");

    Info(InfoChar0GaussLinearEq, 5,
         "MakeIntSystem2: lcm: ", lcm);

    system.lcm := lcm;
    system.int_mat := system.mat * lcm;
    system.int_vecs := system.vecs * lcm;

    system.p := p;
    system.precision := precision;
    system.padic_family := PurePadicNumberFamily(p, precision);

    system.mat_mod_p := system.int_mat * Z(system.p)^0;
    system.mat_mtx64 := MTX64_Matrix(system.mat_mod_p);
    # TODO: Benchmark impact of this
    # ConvertToMatrixRep(system.mat_mod_p);

    system.mtx64 := MTX64_Echelize(system.mat_mtx64);

    # these are the variables we're interested in solving for
    n := MTX64_Matrix_NumRows(system.mtx64.remnant);
    system.solvable_variables := MTX64_EmptyBitString(n);
    for r in [0..n-1] do
        if MTX64_DNzl(system.mtx64.remnant, r) = fail then
            MTX64_SetEntryOfBitString(system.solvable_variables, r);
        fi;
    od;
    system.solvable_variable_poss := MTX64_PositionsBitString(system.solvable_variables);
    system.nr_solvable_variables := Length(system.solvable_variable_poss);

    # Some things we have to do manually atm
    system.col_select := MTX64_PositionsBitString(system.mtx64.colSelect);
    system.row_select := MTX64_PositionsBitString(system.mtx64.rowSelect);

    return system;
end);

InstallGlobalFunction( C0GAUSS_SolutionIntMatVecs_Padic,
function(system)
    local p
        , mat_int
        , sol_p
        , sol_int
        , sol_padic
        , residue
        , residue_p
        , mult
        , ppower
        , iterations
        , done;

    iterations := 0;
    ppower := 1;
    done := false;

    p := system.p;

    mult := system.mtx64.multiplier;
    mat_int := system.int_mat{ MTX64_PositionsBitString(system.mtx64.rowSelect) }{ MTX64_PositionsBitString(system.mtx64.colSelect) };
    residue := system.int_vecs{ MTX64_PositionsBitString(system.mtx64.rowSelect) };
    sol_padic := [];

    while not IsZero(residue) and
        (iterations < system.precision) do
        residue_p := MTX64_Matrix(residue * Z(p)^0);
        sol_p := mult * residue_p;
        sol_int := MatIntFFESymm(MTX64_ExtractMatrix(sol_p));
        sol_padic := sol_padic - ppower * sol_int{system.solvable_variable_poss};

        residue := (residue + mat_int * sol_int) / p;
        iterations := iterations + 1;
        ppower := p * ppower;
    od;
    system.sol_padic := sol_padic;
    system.sol_rat := MatRationalReconstruction(ppower, sol_padic);
end);

InstallGlobalFunction( C0GAUSS_SolutionMat_Padic,
function(mat, vecs, opts...)
    local system, prime, precision;

    prime := C0GAUSS_Padic_Prime;
    precision := C0GAUSS_Padic_Precision;

    if Length(opts) = 1 then
        if IsBound(opts[1].prime) then
            prime := opts[1].prime;
        fi;
        if IsBound(opts[1].precision) then
            precision := opts[1].precision;
        fi;
    fi;

    system := C0GAUSS_SetupMatVecsSystem_Padic( mat, vecs, prime, precision);
    C0GAUSS_SolutionIntMatVecs_Padic(system);

    return system.sol_rat;
end);

InstallGlobalFunction( C0GAUSS_SolutionMatVecs_Padic,
function(mat, vec)
    local vi, system, res, t, t2, tmpmat, tmpvecs;

    res := rec();
    res.solutions := [];
    res.mat := mat;
    res.vec := vec;

    t := NanosecondsSinceEpoch();
    system := C0GAUSS_SetupMatVecsSystem_Padic( TransposedMatMutable(mat), TransposedMat(vec)
                                                 , C0GAUSS_Padic_Prime
                                                 , C0GAUSS_Padic_Precision );
    t :=  NanosecondsSinceEpoch() - t;
    Info(InfoChar0GaussLinearEq, 1, "setup took: ", t/1000000., " msec\n");
    if Length(system.solvable_variables) > 0 then
        Info(InfoChar0GaussLinearEq, 5,
             "Solving for: ", Length(system.transposed_vecs), " rhs\n");
        t2 := NanosecondsSinceEpoch();
        for vi in [1..Length(system.transposed_vecs)] do
            t := NanosecondsSinceEpoch();
            C0GAUSS_SolutionIntMatVec_Padic(system, vi);
            t := NanosecondsSinceEpoch() - t;
            Info(InfoChar0GaussLinearEq, 1, "solving rhs took: ", t/1000000., " msec\n");

            Add(res.solutions, system.rat_solution);
        od;
        t2 := NanosecondsSinceEpoch() - t2;
        Info(InfoChar0GaussLinearEq, 1, "solving all rhs took: ", t2/1000000., " msec\n");


        res.solutions := TransposedMatMutable(res.solutions);
    fi;

    # Debugging
    # res.system := system;
    return res;
end);
