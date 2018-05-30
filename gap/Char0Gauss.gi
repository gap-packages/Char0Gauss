#
# Char0Gauss: Linear algebra stuff
#
# Implementations
#
InstallGlobalFunction( Char0Gauss_Example,
function()
	Print( "This is a placeholder function, replace it with your own code.\n" );
end );

#
# Exists as "RatNumberFromModular" in the edim package as well.
#
InstallGlobalFunction( FareyReconstruction,
function(N, r)
    local a0, b0, a1, b1, a2, b2, q;

    a0 := N;
    b0 := 0;
    a1 := r;
    b1 := 1;

    while 2 * a1^2 > N - 1 do
        q := QuoInt(a0, a1);

        a2 := a0 - q * a1;
        b2 := b0 - q * b1;

        a0 := a1; b0 := b1;
        a1 := a2; b1 := b2;
    od;
    if (2 * b1^2 <= N - 1) and GcdInt(a1, b1) = 1 then
        return a1 / b1;
    fi;

    return fail;
end );

