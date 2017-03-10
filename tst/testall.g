#
# Char0Gauss: Linear algebra stuff
#
# This file runs package tests. It is also referenced in the package
# metadata in PackageInfo.g.
#
LoadPackage( "Char0Gauss" );

TestDirectory(DirectoriesPackageLibrary( "Char0Gauss", "tst" ),
  rec(exitGAP := true));

FORCE_QUIT_GAP(1); # if we ever get here, there was an error
