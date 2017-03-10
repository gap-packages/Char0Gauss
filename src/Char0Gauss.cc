/*
 * Char0Gauss: Linear algebra stuff
 */

// Include gmp.h *before* switching to C mode, because GMP detects when
// compiled from C++ and then does some things differently, which would
// cause an error if called from within extern "C". But libsing.h
// (indirectly) includes gmp.h ...
#include <gmp.h>

extern "C" {
#include "src/compiled.h"          /* GAP headers                */
    #include "src/rational.h"

    Int EqRat(Obj l, Obj r);
    Obj ZeroRat(Obj x);
    Obj AInvRat(Obj x);
    Obj InvRat(Obj x);
    Obj ProdRat(Obj l, Obj r);
    Obj SumRat(Obj l, Obj r);
    Obj DiffRat(Obj l, Obj r); // Hack, because AInvRat does not accept integers..
}

/*
 * Fairly naive (but hopefully faster) gaussian elimination for
 * matrices of rationals
 * Input is assumed to be a list of lists of correct shape, not empty, containing only rationals,
 * and all the rest breakage of it
 *
 * This function will destroy its input and if you care about it save
 * it ahead of time!
 */

// dst += src * a
static inline Obj AddRationalRows(Obj dst, Obj src, Obj a)
{
    Obj tmp;
    UInt i;

    for(i=1; i <= LEN_PLIST(dst); i++) {
        tmp = ELM_PLIST(src, i);
        tmp = ProdRat(tmp, a);
        tmp = SumRat(ELM_PLIST(dst,i),tmp);
        SET_ELM_PLIST(dst, i, tmp);
        CHANGED_BAG(dst);
    }
    return 0L;
}

static inline Obj MultRationalRow(Obj dst, Obj a)
{
    Obj tmp;
    UInt i;

    for(i=1;i<=LEN_PLIST(dst);i++) {
        tmp = ProdRat(ELM_PLIST(dst, i), a);
        SET_ELM_PLIST(dst, i, tmp);
        CHANGED_BAG(dst);
    }
    return 0L;
}

Obj SemiEchelonRationals(Obj self, Obj mat)
{
    Int p;
    UInt i, j;
    UInt nrows, ncols;
    UInt nvecs, nheads, nnzheads;

    Obj row;
    Obj vectors;
    Obj heads;
    Obj result;
    Obj nzheads;
    Obj x;
    Obj zero;
    Obj inv;
    Obj tmp, tmp2;

    // ZeroRat ignores its argument...
    zero = ZeroRat(inv);

    nrows = LEN_PLIST(mat);
    ncols = LEN_PLIST(ELM_PLIST(mat, 1));

    heads = NEW_PLIST(T_PLIST, ncols);
    nheads = 1;
    SET_LEN_PLIST(heads, 0);

    vectors = NEW_PLIST(T_PLIST, ncols);
    nvecs = 1;
    SET_LEN_PLIST(vectors, 0);

    nzheads = NEW_PLIST(T_PLIST, ncols);
    nnzheads = 1;
    SET_LEN_PLIST(nzheads, 0);

    for(i=1;i<=nrows;i++) {
        row = ELM_PLIST(mat, i);

        for(j=1;j<nnzheads;j++) {
            p = INT_INTOBJ(ELM_PLIST(nzheads, j));

            x = ELM_PLIST(row, p);
            if (EqRat(x, zero) == 0) {
                tmp = ELM_PLIST(vectors, j);
                tmp2 = DiffRat(zero, x);
                AddRationalRows(row, tmp, tmp2);
            }
        }

        /* First non-zero position */
        for(j=1;EqRat(ELM_PLIST(row,j), zero);j++) {
        };
        if(j<=ncols){
            /* New basis vector */

            /* This should, can't and will not fail. It will NOT FAIL */
            inv = InvRat(ELM_PLIST(row,j));
            MultRationalRow(row, inv);
            SET_LEN_PLIST(vectors, nvecs);
            SET_ELM_PLIST(vectors, nvecs++, row);
            CHANGED_BAG(vectors);
            SET_LEN_PLIST(nzheads, nnzheads);
            SET_ELM_PLIST(nzheads, nnzheads++, INTOBJ_INT(j));
            CHANGED_BAG(nzheads);
            SET_LEN_PLIST(heads, j);
            SET_ELM_PLIST(heads, j, INTOBJ_INT(nvecs-1));
            CHANGED_BAG(heads);
        }
    }

    result = NEW_PREC(2);
    AssPRec(result, RNamName("heads"), heads);
    AssPRec(result, RNamName("vectors"), vectors);

    return result;
}

typedef Obj (* GVarFunc)(/*arguments*/);

#define GVAR_FUNC_TABLE_ENTRY(srcfile, name, nparam, params) \
  {#name, nparam, \
   params, \
   (GVarFunc)name, \
   srcfile ":Func" #name }

// Table of functions to export
static StructGVarFunc GVarFuncs [] = {
    GVAR_FUNC_TABLE_ENTRY("Char0Gauss.c", SemiEchelonRationals, 1, "mat"),

	{ 0 } /* Finish with an empty entry */

};

/******************************************************************************
*F  InitKernel( <module> )  . . . . . . . . initialise kernel data structures
*/
static Int InitKernel( StructInitInfo *module )
{
    /* init filters and functions                                          */
    InitHdlrFuncsFromTable( GVarFuncs );

    /* return success                                                      */
    return 0;
}

/******************************************************************************
*F  InitLibrary( <module> ) . . . . . . .  initialise library data structures
*/
static Int InitLibrary( StructInitInfo *module )
{
    /* init filters and functions */
    InitGVarFuncsFromTable( GVarFuncs );

    /* return success                                                      */
    return 0;
}

/******************************************************************************
*F  InitInfopl()  . . . . . . . . . . . . . . . . . table of init functions
*/
static StructInitInfo module = {
 /* type        = */ MODULE_DYNAMIC,
 /* name        = */ "Char0Gauss",
 /* revision_c  = */ 0,
 /* revision_h  = */ 0,
 /* version     = */ 0,
 /* crc         = */ 0,
 /* initKernel  = */ InitKernel,
 /* initLibrary = */ InitLibrary,
 /* checkInit   = */ 0,
 /* preSave     = */ 0,
 /* postSave    = */ 0,
 /* postRestore = */ 0
};

extern "C"
StructInitInfo *Init__Dynamic( void )
{
    return &module;
}
