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
                 , number_variables := Ncols(mat)
                 , number_equations := Nrows(mat) );

    t := NanosecondsSinceEpoch();
    lcm := C0GAUSS_LeastCommonDenominator(mat);
    lcm := LcmInt(lcm, C0GAUSS_LeastCommonDenominator(vecs));
    t := NanosecondsSinceEpoch() - t;
    Info(InfoChar0GaussLinearEq, 1, "computing the LCM took: ", t / 1000000., " msec");

    system.lcm := lcm;
    system.int_mat := system.mat * lcm;
    system.int_vecs := system.vecs * lcm;

    system.p := p;
    system.precision := precision;
    system.padic_family := PurePadicNumberFamily(p, precision);

    t := NanosecondsSinceEpoch();
    system.mtx64 := system.int_mat * Z(p)^0;
    system.mtx64!.ring := GF(p);
    t := NanosecondsSinceEpoch() - t;
    Info(InfoChar0GaussLinearEq, 1, "reduction mod p took: ", t / 1000000., " msec");


    t := NanosecondsSinceEpoch();
    system.mtx64 := ConvertSparseMatrixToMatrix(system.mtx64);
    t := NanosecondsSinceEpoch() - t;
    Info(InfoChar0GaussLinearEq, 1, "sparse matrix to matrix took: ", t / 1000000., " msec");


    t := NanosecondsSinceEpoch();
    system.mtx64 := MTX64_Matrix(system.mtx64);
    t := NanosecondsSinceEpoch() - t;
    Info(InfoChar0GaussLinearEq, 1, "MTX64 creation took: ", t / 1000000., " msec");


    t := NanosecondsSinceEpoch();
    system.mtx64 := MTX64_Echelize(system.mtx64);
    t := NanosecondsSinceEpoch() - t;
    Info(InfoChar0GaussLinearEq, 1, "echelize took: ", t / 1000000., " msec");


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

    if Length(MTX64_PositionsBitString(system.solvable_variables)) = 0 then
        system.sol_padic := [];
        system.sol_rat := [];
        return;
    fi;

    iterations := 0;
    ppower := 1;
    done := false;

    p := system.p;

    mult := system.mtx64.multiplier;
    mat_int := CertainColumns(CertainRows(system.int_mat, system.row_select), system.col_select);
    residue := CertainRows(system.int_vecs, system.row_select);
    sol_padic := SparseZeroMatrix(Length(system.solvable_variable_poss), Ncols(system.int_vecs), Rationals);

    done := false;

    while not (IsSparseZeroMatrix(residue) or done) do
        # iterations < system.precision) do
        residue_p := MTX64_Matrix(ConvertSparseMatrixToMatrix(residue) * Z(p)^0);
        sol_p := mult * residue_p;
        sol_int := SparseMatrix(MatIntFFESymm(MTX64_ExtractMatrix(sol_p)), Rationals);
        sol_padic := sol_padic - ppower * CertainRows(sol_int, system.solvable_variable_poss);

        residue := (residue + mat_int * sol_int) * (1 / p);
        iterations := iterations + 1;
        ppower := p * ppower;

        if (iterations > system.precision) and (iterations mod 100 = 0) then
            Info(InfoChar0GaussLinearEq, 10, iterations, " iterations done.\n");
            system.sol_rat := MatRationalReconstruction(ppower, ConvertSparseMatrixToMatrix(sol_padic));
            if PositionProperty(Concatenation(system.sol_rat), x -> x = fail) = fail then
                done := true;
            fi;
        fi;
    od;
    system.sol_padic := sol_padic;
#     system.sol_rat := MatRationalReconstruction(ppower, sol_padic);

    if PositionProperty(Concatenation(system.sol_rat), x -> x = fail) <> fail then
        Error("failed to recosntruct");
    fi;
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
function(res)
    local vi, system, t, t2, tmpmat, tmpvecs, sv;

    res.solutions := [];

    t := NanosecondsSinceEpoch();
    system := C0GAUSS_SetupMatVecsSystem_Padic( res.mat, res.vec
                                                 , C0GAUSS_Padic_Prime
                                                 , C0GAUSS_Padic_Precision );
    t := NanosecondsSinceEpoch() - t;
    Info(InfoChar0GaussLinearEq, 1, "setup took: ", t/1000000., " msec");

    t := NanosecondsSinceEpoch();
    C0GAUSS_SolutionIntMatVecs_Padic(system);
    t := NanosecondsSinceEpoch() - t;
    Info(InfoChar0GaussLinearEq, 1, "solving took: ", t/1000000., " msec");
    Info(InfoChar0GaussLinearEq, 1, "number of solvables: ", Length(system.solvable_variables), ".");

    # TODO: Document why this is correct
    sv := system.col_select{ MTX64_PositionsBitString(system.solvable_variables) };

    res.solutions := ListWithIdenticalEntries( Ncols(res.mat), fail );
    res.solutions{sv} := List(system.sol_rat, x->SparseMatrix([x],Rationals));

    res.mat := CertainRows(system.int_mat, system.row_select);
    res.vec := CertainRows(system.int_vecs, system.row_select);

    return res;
end);
