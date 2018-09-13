#
gap> mat := [[1/3,1/4],[3,2]];;
gap> C0GAUSS_LeastCommonDenominator(mat);
12
gap> C0GAUSS_LeastCommonDenominator(mat, 5);
60
gap> mats := SparseMatrix(mat, Rationals);;
gap> C0GAUSS_LeastCommonDenominator(mats);
12
gap> C0GAUSS_LeastCommonDenominator(mats, 5);
60

#
