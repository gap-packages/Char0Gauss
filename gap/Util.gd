DeclareGlobalFunction( "C0GAUSS_FoldList2" );
DeclareGlobalFunction( "MatIntFFESymm" );

DeclareOperation( "C0GAUSS_LeastCommonDenominator", [IsObject]); # There is no more specialised filter that captures
# sparse matrices and matrices...
DeclareOperation( "C0GAUSS_LeastCommonDenominator", [IsObject, IsPosInt]);


