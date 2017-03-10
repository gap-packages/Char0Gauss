#
# Char0Gauss: Linear algebra stuff
#
# Reading the declaration part of the package.
#
_PATH_SO:=Filename(DirectoriesPackagePrograms("Char0Gauss"), "Char0Gauss.so");
if _PATH_SO <> fail then
    LoadDynamicModule(_PATH_SO);
fi;
Unbind(_PATH_SO);

ReadPackage( "Char0Gauss", "gap/Char0Gauss.gd");
