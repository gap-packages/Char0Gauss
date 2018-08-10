#
# Solving linear equations over the integers/rationals by Dixon/Hensel lifting
#
DeclareInfoClass("InfoChar0GaussPadics");
DeclareInfoClass("InfoChar0GaussLinearEq");

# Setup status structure
DeclareGlobalFunction( "C0GAUSS_SetupMatVecsSystem_Padic" );

DeclareGlobalFunction( "C0GAUSS_Padic_Presolve" );

# TODO: Should this go into MeatAxe64?
DeclareGlobalFunction( "MTX64_SolutionsMatWithEchelized" );
DeclareGlobalFunction( "MTX64_SolutionsMatWithEchelizedSelected" );

# solve for one vector (integer version)
DeclareGlobalFunction( "C0GAUSS_SolutionIntMatVecs_Padic" );

# solve for a list of vectors (rationals)
DeclareGlobalFunction( "MAJORANA_SolutionMatVecs_Padic" );
DeclareGlobalFunction( "C0GAUSS_SolutionMatVec_Padic" );

# Some default values, used in MAJORANA_SolutionMatVecs_Padic
BindGlobal( "MAJORANA_Padic_Prime", 191 );
BindGlobal( "MAJORANA_Padic_Precision", 100 );
BindGlobal( "MAJORANA_Padic_Iterations", 100 );

# Some default values
BindGlobal( "C0GAUSS_Padic_Prime", 191 );
BindGlobal( "C0GAUSS_Padic_Precision", 100 );
BindGlobal( "C0GAUSS_Padic_Iterations", 100 );

