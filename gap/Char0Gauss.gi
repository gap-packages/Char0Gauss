#
# Char0Gauss: Linear algebra stuff
#
# Implementations
#

LoadPackage("meataxe64");

#
# Exists as "RatNumberFromModular" in the edim package as well.
#
RationalReconstruction := function(N, r)
    local a0, b0, a1, b1, a2, b2, q;

    a0 := N;
    b0 := 0;
    a1 := r;
    b1 := 1;

    while 2 * a1^2 > N - 1 do
        q := QuoInt(a0, a1);

        a2 := a0 - q * a1;
        b2 := b0 - q * b1;

        a0 := a1; b0 := b1;
        a1 := a2; b1 := b2;
    od;
    if (2 * b1^2 <= N - 1) and GcdInt(a1, b1) = 1 then
        return a1 / b1;
    fi;

    return fail;
end;

# testing the thing
Display(RationalReconstruction(101, 77));
for n in [10..1000] do
    count := 0;
    for i in [0..n-1] do
        q := RationalReconstruction(n, i);
        if q = fail then
        else
    #        Print(count, " : ", i, ": ", q, "\n");
            count := count + 1;
        fi;
    od;
    ratio := 1. * (count / n);
    Print("[0..", n, ") size: ", count, " ratio: ", ratio, "\n");
od;

#ErrorTolerantLifting := function(N, r)
#    local a0, b0, a1, b1, q;
#    a0 := N;
#    b0 := 0;
#    a1 := r;
#    b1 := 1;
#    repeat
#        q := [
#            a0 - Floor(LCM_INT(a0, b0) / a1^2) * a1,
#            b0 - Floor(LCM_INT(a1, b1) / b1^2) * b1
#        ];
#        a0 := a1;
#        b0 := b1;
#        a1 := q[1];
#        b1 := q[2];
#    until a1^2 + b1^2 >= a0^2 + b0^2;
#    if a0^2 + b0^2 < N then
#        return a0/b0;
#    fi;
#    return false;
#end;

Char0Gauss := function(A, b)
    local s, solns, primes, v, p, i, j, t;
    # scan the whole matrix to choose the primes
    primes := [101, 239, 997, 5813, 37199];
    for i in [1..Length(primes)] do
    od;
    # solution vector
    solns := 0 * [1..Length(primes)];
    for i in [1..Length(primes)] do
        p := primes[i];
        Display(p);
        solns[i] := SolutionMat(One(Z(p)) * A, One(Z(p)) * b);
    od;
    Print("solns: ", solns, "\n");
    # dirty and incorrect way to get the numbers back
    v := 0 * [1..Length(solns[1])];
    # iterate vectors corresponding to other primes
    for i in [1..Length(solns[1])] do
      v[i] := 0 * [1..Length(solns)];
      for j in [1..Length(solns)] do
        v[i][j] := Int(solns[j][i]);
      od;
      v[i] := RationalReconstruction(Product(primes), ChineseRem(primes, v[i]));
    od;
    return v;
end;

A := [
  [1/2, 2, 3],
  [0, 2, 5],
  [0, 0, 4576/18734]
];
b := [1, 1, 1];
Display(A);
Display(SolutionMat(A, b));
Display(Char0Gauss(A, b));
