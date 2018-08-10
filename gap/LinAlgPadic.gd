#
# Solving linear equations over the integers/rationals by Dixon/Hensel lifting
#
DeclareInfoClass("InfoChar0GaussPadics");
DeclareInfoClass("InfoChar0GaussLinearEq");

# Setup status structure
DeclareGlobalFunction( "C0GAUSS_SetupMatVecsSystem_Padic" );

DeclareGlobalFunction( "C0GAUSS_SolutionIntMatVecs_Padic" );
DeclareGlobalFunction( "C0GAUSS_SolutionIntMatVec_Padic" );

# solve for a list of vectors (rationals)
DeclareGlobalFunction( "C0GAUSS_SolutionMatVecs_Padic" );
DeclareGlobalFunction( "C0GAUSS_SolutionMatVec_Padic" );
DeclareGlobalFunction( "C0GAUSS_SolutionMat_Padic" );

# Some default values
BindGlobal( "C0GAUSS_Padic_Prime", 191 );
BindGlobal( "C0GAUSS_Padic_Precision", 100 );
