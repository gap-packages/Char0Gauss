gap> genmat_random := function(m, n)
>     return [
>         RandomMat(m, n, Rationals),
>         RandomMat(1, n, Rationals)[1]
>     ];
> end;;
gap> genmat_invertible := function(m, n)
>     local A, b;
>     if m * n < 2^12 then
>         if m = n then
>             A := RandomInvertibleMat(n, Rationals);
>         else
>             repeat
>                 A := RandomMat(m, n, Rationals);
>             until Rank(A) = Minimum(m, n);
>         fi;
>     else
>         A := RandomMat(m, n, Rationals);
>     fi;
>     b := RandomMat(1, n, Rationals)[1];
>     return [A, b];
> end;;
gap> genmat_invertible_small_sln := function(m, n)
>     local A, v, bound;
>     if m * n < 2^12 then
>         if m = n then
>             A := RandomInvertibleMat(n, Rationals);
>         else
>             repeat
>                 A := RandomMat(m, n, Rationals);
>             until Rank(A) = Minimum(m, n);
>         fi;
>     else
>         A := RandomMat(m, n, Rationals);
>     fi;
>     bound := (m*n)^2;
>     v := List([1..n], x -> Random([-bound..bound]) / Random([1..bound]));
>     return [A, v * A];
> end;;
gap> run_tests_char0gauss := function(m, n, gen_mat, times)
>     local obj, A, b, v, w, u, i, times1, times2, times3, t1, t2, t3, small_primes;
>     times1 := [];
>     times2 := [];
>     times3 := [];
>     for i in [1..times] do
>         obj := gen_mat(m, n);
>         A := obj[1];
>         b := obj[2];
>         t2 := NanosecondsSinceEpoch();
>         w := Char0Gauss_SolutionMat(A, b, [], false);
>         t2 := NanosecondsSinceEpoch() - t2;
>         t1 := NanosecondsSinceEpoch();
>         v := SolutionMat(A, b);
>         t1 := NanosecondsSinceEpoch() - t1;
>         Append(times1, [t1 / 1000000000.]);
>         Append(times2, [t2 / 1000000000.]);
>         if not v = w then
>             Display("fail");
>             Display(v);
>             Display(w);
>             return false;
>         fi;
>     od;
>     return true;
> end;;

#
gap> SetInfoLevel(InfoChar0GaussChineseRem, 1);
gap> for i in [1..10] do
>     if run_tests_char0gauss(i, i, genmat_invertible, 10) = fail then
>         break;
>     fi;
> od;;
