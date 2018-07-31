#
# Solving linear equations over the integers/rationals by Dixon/Hensel lifting
#
DeclareInfoClass("InfoMajoranaPadics");
DeclareInfoClass("InfoMajoranaLinearEq");

# Setup status structure
DeclareGlobalFunction( "MAJORANA_SetupMatVecsSystem_Padic" );

DeclareGlobalFunction( "MAJORANA_Padic_Presolve" );

# solve for one vector (integer version)
DeclareGlobalFunction( "MAJORANA_SolutionIntMatVec_Padic" );

# solve for a list of vectors (rationals)
DeclareGlobalFunction( "MAJORANA_SolutionMatVecs_Padic" );
DeclareGlobalFunction( "C0GAUSS_SolutionMatVec_Padic" );

# Some default values, used in MAJORANA_SolutionMatVecs_Padic
BindGlobal( "MAJORANA_Padic_Prime", 191 );
BindGlobal( "MAJORANA_Padic_Precision", 100 );
BindGlobal( "MAJORANA_Padic_Iterations", 100 );

# FIXME: Todo?
# DeclareGlobalFunction( "MAJORANA_NullspaceIntMat_Padic" );
# DeclareGlobalFunction( "MAJORANA_NullspaceMat_Padic" );

