# Solving linear equations over the integers/rationals by Dixon/Hensel lifting
#
# This is a GAP prototype which already works quite a bit faster than any code
# inside GAP, for reasonably large systems (thousands of equations)
#
# If you find any bugs, email markus.pfeiffer@morphism.de
#
# TODO:
# * Carry denominator forward
# * actually only solve the solvable variables.
# * Make a better implementation of the padics code. Its currently pretty brittle
#   and hacky
# * More tests
# * look at flint that has some of this functionality
# * Implement in C or Rust (or Julia)?
# * Parallelisation strategies
# * use meataxe64
#

# Just to make sure we're not shooting ourselves
# in the foot with inconsistent entries.
CheckSystem := function(system)
    local b, r, c;

    Info( InfoMajoranaLinearEq, 5,
          " testing system of equation structure" );
    if not IsPrime(system.p) then
        Error("p is not prime");
    fi;

    for r in system.int_mat do
        for c in r do
            if DenominatorRat(c) <> 1 then
                Error("Non-1-denominator in system.int_mat");
            fi;
        od;
    od;

    for r in system.int_vecs do
        for c in r do
            if DenominatorRat(c) <> 1 then
                Error("Non-1-denominator in system.int_vecs");
            fi;
        od;
    od;
    Info( InfoMajoranaLinearEq, 5,
          " success.");
end;

InstallGlobalFunction( MAJORANA_Padic_Presolve,
function(system)
    local v, r, i;

    system.mat_mod_p := system.int_mat * Z(system.p)^0;
    ConvertToMatrixRep(system.mat_mod_p);

    system.echelon := EchelonMatTransformation(system.mat_mod_p);

    # FIXME: refactor
    v := List( system.echelon.vectors
             , x -> PositionsProperty(x, y -> not IsZero(y) ) );
    system.solvable_rows := PositionsProperty( v, z -> Length(z) = 1);

    system.interesting_rows := Set(Concatenation( List( system.echelon.coeffs,
                                               r -> PositionsProperty( r, x -> not IsZero(x)) ) ));
    system.solvable_variables := Concatenation( Filtered( v
                                                        , z -> Length(z) = 1) );
    system.unsolvable_variables := Difference([1..system.number_variables], system.solvable_variables);
    system.unsolvable_rows := Difference([1..system.number_equations], system.solvable_rows);
    system.uninteresting_rows := Difference([1..system.number_equations], system.interesting_rows);

    return system;
end );

InstallGlobalFunction(MAJORANA_SetupMatVecsSystem_Padic,
function(mat, vecs, p, precision, max_iter)
    local system, mmults, vmults, lcm;
    system := rec( mat := mat
                 , vecs := vecs
                 , number_variables := Length(mat[1])
                 , number_equations := Length(mat) );

    #  MakeIntSystem(system);
    Info(InfoMajoranaLinearEq, 5,
         "MakeIntSystem2: computing denominator lcms" );

    mmults := List(system.mat, x -> C0GAUSS_FoldList2(x, DenominatorRat, LcmInt));
    vmults := List(system.vecs, x -> C0GAUSS_FoldList2(x, DenominatorRat, LcmInt));
    lcm := C0GAUSS_FoldList2(Concatenation(mmults, vmults), IdFunc, LcmInt);

    Info(InfoMajoranaLinearEq, 5,
         "MakeIntSystem2: lcm: ", lcm);

    system.lcm := lcm;
    system.int_mat := system.mat * lcm;
    system.int_vecs := system.vecs * lcm;

    system.p := p;
    system.precision := precision;
    system.padic_family := PurePadicNumberFamily(p, precision);
    system.padic_iterations := max_iter;

    MAJORANA_Padic_Presolve(system);

    # Transposingpalooza
    # FIXME: cleanup and only keep the ones we need
    system.transposed_int_mat := TransposedMat(system.int_mat);
    system.transposed_coeffs := TransposedMat(system.echelon.coeffs);
    system.transposed_vecs := TransposedMat(system.int_vecs);
    system.lifted_coeffs := List(system.transposed_coeffs, y -> List(y, IntFFESymm));

    system.solution_denominator := 1;

    return system;
end);

InstallGlobalFunction( MAJORANA_SolutionIntMatVec_Padic,
function(system, vi)
    local
        p,

        # These are *integer* vectors
        tmp_soln, soln,
        residue,

        # These are vectors in GF(p)
        vec_p,
        soln_p, soln_pp,

        done, iterations,
        soln_padic,
        ppower, sol, x, y, i, t, t2,
        k, old_denom, denom, vecd, iter,
        recovered_rat, recovered_rat2, ll, zero_p;

    p := system.p;

    # Accumulator for integer solution
    # FIXME: only solved variables? Should cut down on memory use
    #        and how many p-adic expansions we have to actually turn
    #        into denominators
    soln := ListWithIdenticalEntries(system.number_variables, 0);

    # These are the *integer* residuals of the RHS
    # initially this is the RHS we're solving for
    # FIXME: for testing, only do one RHS right now
    #        this will trivially generalise fortunately
    #        but there might be a point in not solving all
    #        RHS at the same time, in case we discover enough 
    #        of the denominator to not have to approximate?
    residue := MutableCopyMat(system.transposed_vecs[vi])
               * system.solution_denominator;

    done := false;
    iterations := 0;
    ppower := 1; # this is p^iterations, FIXME: can we avoid this?

    # digits in the p-adic approximation to the solution, stored as
    # an integer
    soln_padic := ListWithIdenticalEntries(system.number_variables, 0);
    soln_p := ListWithIdenticalEntries(system.number_variables, Zero(Z(p)));

    t2 := NanosecondsSinceEpoch();
    while (not done) do
        iterations := iterations + 1;

        if iterations mod 100 = 0 then
            Info(InfoMajoranaLinearEq, 10, STRINGIFY(iterations, " iterations"));
        fi;

        # solve the system mod p
        vec_p := Z(p)^0 * residue;
        soln_p := vec_p * system.transposed_coeffs;
        soln_p := List( system.echelon.heads,
                      function(x)
                          if x > 0 then
                              return soln_p[x];
                          else
                              return Zero(soln_p[1]);
                          fi;
                      end);

        # Convert the solution from GF(p) to integers -p/2..p/2-1
        y := List(soln_p, IntFFESymm);
        AddRowVector(soln, y, ppower);

        y := -y * system.transposed_int_mat;
        
        for i in system.interesting_rows do
            residue[i] := (residue[i] + y[i])/p;
        od;

       # residue{ system.interesting_rows } := residue{ system.interesting_rows } / p;

        ppower := ppower * p;

        Info(InfoMajoranaLinearEq, 10, "soln:    ", soln);
        Info(InfoMajoranaLinearEq, 10, "y:       ", y);
        Info(InfoMajoranaLinearEq, 10, "residue: ", residue);

        # Solution found?
        if IsZero( residue ) then
            Info(InfoMajoranaLinearEq, 1,
                 "found an integer solution");

            # FIXME: I don't like this state struct design at the moment
            system.int_solution := soln;
            system.rat_solution := soln;
            return true;
        else
            if iterations >= system.precision then
                t2 := NanosecondsSinceEpoch() - t2;
                Info(InfoMajoranaLinearEq, 1, "\n",
                     "reaching iteration limit", iterations, " took, ", t2 / 1000000., " msec");
                Info(InfoMajoranaLinearEq, 1,
                     "average bits          ", Length(soln), " ", Float(Sum(List(soln, Log2Int)) / Length(soln)), "\n",
                     "average bits (sol)    ", Length(system.solvable_variables), " ", Float(Sum(List(soln{system.solvable_variables}, Log2Int)) / (Length(system.solvable_variables) + 1)), "\n",
                     "average bits (unsol)  ", Length(system.unsolvable_variables), " ", Float(Sum(List(soln{system.unsolvable_variables}, Log2Int)) / (Length(system.unsolvable_variables) + 1)), "\n",
                     "average bits (resid)  ", Length(residue), " ", Float(Sum(List(soln, Log2Int)) / Length(soln)), "\n",
                     "average bits (intr)   ", Length(system.interesting_rows), " ", Float(Sum(List(residue{system.interesting_rows}, Log2Int)) / (Length(system.interesting_rows) + 1)), "\n",
                     "average bits (unintr) ", Length(system.uninteresting_rows), " ", Float(Sum(List(residue{system.uninteresting_rows}, Log2Int)) / (Length(system.uninteresting_rows) + 1)), "\n"
                    );

                Info(InfoMajoranaLinearEq, 5,
                     "reached iteration limit, trying to recover rational solution");

                # Recover rational solutions
                recovered_rat := List( soln{ system.solvable_variables }
                                     , x -> RationalReconstruction(ppower, x) );
                if ForAny(recovered_rat, x -> x = fail) then
                    Error("failed to recover");
                else
                    system.rat_solution := ListWithIdenticalEntries( system.number_variables, 0 );
                    system.rat_solution{ system.solvable_variables } := recovered_rat;
                    return true;
                fi;
            fi;
        fi;
    od;
end );

# FIXME: This has to be compatible with MAJORANA_SolutionMatVecs
# FIXME: Make compatible with sparse matrices.
InstallGlobalFunction( MAJORANA_SolutionMatVecs_Padic,
function(mat, vecs)
    local vi, system, res, t, t2, tmpmat, tmpvecs;

    res := rec();
    res.solutions := [];
    res.mat := mat;
    res.vec := vecs;

    mat := ConvertSparseMatrixToMatrix(mat);
    vecs := ConvertSparseMatrixToMatrix(vecs);

    t := NanosecondsSinceEpoch();
    system := MAJORANA_SetupMatVecsSystem_Padic( MutableCopyMat(mat), vecs
                                                 , MAJORANA_Padic_Prime
                                                 , MAJORANA_Padic_Precision
                                                 , MAJORANA_Padic_Iterations );
    t :=  NanosecondsSinceEpoch() - t;
    Info(InfoMajoranaLinearEq, 1, "setup took: ", t/1000000., " msec\n");
    if Length(system.solvable_variables) > 0 then
        Info(InfoMajoranaLinearEq, 5,
             "Solving for: ", Length(system.transposed_vecs), " rhs\n");
        t2 := NanosecondsSinceEpoch();
        for vi in [1..Length(system.transposed_vecs)] do
            t := NanosecondsSinceEpoch();
            MAJORANA_SolutionIntMatVec_Padic(system, vi);
            t := NanosecondsSinceEpoch() - t;
            Info(InfoMajoranaLinearEq, 1, "solving rhs took: ", t/1000000., " msec\n");

            Add(res.solutions, system.rat_solution);
        od;
        t2 := NanosecondsSinceEpoch() - t2;
        Info(InfoMajoranaLinearEq, 1, "solving all rhs took: ", t2/1000000., " msec\n");


        res.solutions := TransposedMatMutable(res.solutions);
        res.solutions := List(res.solutions, x -> SparseMatrix([x], Rationals));
    fi;
    res.solutions{ system.unsolvable_variables } := ListWithIdenticalEntries(Length(system.unsolvable_variables), fail);

#    tmpmat := system.lifted_coeffs * system.mat;
#    tmpvecs := system.lifted_coeffs * system.vecs;

#    res.mat := SparseMatrix(tmpmat{ system.unsolvable_rows }, Rationals);
#    res.vec := SparseMatrix(tmpvecs{ system.unsolvable_rows }, Rationals);

#     res.mat := CertainRows(SparseMatrix(system.mat), system.unsolvable_rows);
#    res.vec := CertainRows(SparseMatrix(system.vecs), system.unsolvable_rows);

    res.mat := SparseMatrix(system.mat{ system.unsolvable_rows }, Rationals );
    res.vec := SparseMatrix(system.vecs{ system.unsolvable_rows }, Rationals );


    res.mat!.ring := Rationals;
    res.vec!.ring := Rationals;

    # Debugging
    res.system := system;
    return res;
end);
