#
gap> testLinAlg := function(n)
>   local A, b, Ac, bc
>        , sol_a, sol_b;
>   A := RandomMat(n,n,Rationals);
>   b := RandomMat(n,1,Rationals);
>   Ac := MutableCopyMat(A);
>   bc := MutableCopyMat(b);
>   sol_a := C0GAUSS_SolutionMat_Padic(Ac, bc);
>   Ac := TransposedMat(MutableCopyMat(A));
>   bc := TransposedMat(MutableCopyMat(b));
>   sol_b := SolutionMat(Ac, bc[1]);
>   return TransposedMat(sol_a)[1] - sol_b;
> end;;

#
gap> IsZero(testLinAlg(30));
true
gap> IsZero(testLinAlg(40));
true
gap> IsZero(testLinAlg(20));
true
gap> IsZero(testLinAlg(25));
true

#
