gap> testLinAlg := function(mat, vecs)
>  local s1, s2;
>  s1 := MAJORANA_SolutionMatVecs(mat, vecs);
>  s2 := MAJORANA_SolutionMatVecs_Padic(mat, vecs);
>  return s1.solutions = s2.solutions;
> end;;

gap> testLinAlg2 := function(mat, vecs)
>  local s1, sls, v, tm, tv;
>  s1 := MAJORANA_SolutionMatVecs_Padic(mat, vecs);
>  sls := [];
>  
>  tm := ConvertSparseMatrixToMatrix(TransposedMat(mat));
>  tv := ConvertSparseMatrixToMatrix(TransposedMat(vecs));
>
>  for v in tv do
>     Add(sls, SolutionMat(tm, v));
>  od;
>
>  return sls;
> end;;




# 
gap> M := SparseMatrix([[1/3, 1],[0,1/2]], Rationals);; b := SparseMatrix([[2],[3]], Rationals);;
gap> testLinAlg(M, b);;

#
gap> M2 := SparseMatrix([[1,1,1,1,1],[0,1,0,1,0],[0,0,1,0,1]], Rationals);; b2 := SparseMatrix([[3],[7],[2]], Rationals);;
gap> testLinAlg(M, b);;


