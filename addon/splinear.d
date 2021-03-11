/**
 * Basic vector, matrix, idx tools.
 *
 * Copyright:   Copyright (C) 2021 by Yinchieh Lai
 * Author:     Yinchieh Lai
 * License:     $(LINK2 http://www.boost.org/LICENSE_1_0.txt, Boost License 1.0)
 */

module addon.splinear;
import addon.tool;
import addon.umfpack;
import addon.spm;

pragma(lib, "umfpack");
//  version(Posix) pragma(link, "amd");
//  pragma(link, "blas");

double[] cr_copyto_ri(Cdouble[] cr, double[] ri = null)
{
    size_t n = cr.length;
    if (ri is null)
        ri.length = n;
    foreach (i, s; cr)
    {
        ri[i] = s.re;
        ri[i + n] = s.im;
    }
    return ri;
}

Cdouble[] ri_copyto_cr(double[] ri, Cdouble[] cr)
{
    size_t n = ri.length / 2;
    if (cr is null)
        cr.length = n;
    size_t j = 0, k = 1;
    double[] cvout = (cast(double*) cr.ptr)[0 .. 2 * n];
    for (size_t i = 0; i < n; ++i)
    {
        cvout[j] = ri[i];
        cvout[k] = ri[i + n];
        j += 2;
        k += 2;
    }
    return cr;
}

T[] linearsolve(T)(Spm!(T) m, T[] b, T[] x = null)
{
    static if (is(T : double))
        return dlinearsolve(m.compact, b, x);
    else static if (is(T : Cdouble))
    {
        if (x is null)
            x.length = b.length;
        double[] bri;
        double[] xri;
        size_t n = b.length;
        bri.length = xri.length = 2 * n;
        cr_copyto_ri(b, bri);
        clinearsolve(m.compact_ri, bri, xri);
        ri_copyto_cr(xri, x);
        bri.destroy();
        xri.destroy();
        return x;
    }
    else
        assert(false, "wrong matrix type");
}

double[] dlinearsolve(SpmCompact!(double) m, double[] bx, double[] x = null)
{
    int n = to!int(bx.length);
    void* Symbolic, Numeric;
    if (x is null)
        x = newvec(n, double.init);
    int[] Ap = m.count;
    int[] Ai = m.id;
    double[] Ax = m.data;
    umfpack_di_symbolic(n, n, Ap.ptr, Ai.ptr, Ax.ptr, &Symbolic, null, null);
    umfpack_di_numeric(Ap.ptr, Ai.ptr, Ax.ptr, Symbolic, &Numeric, null, null);
    umfpack_di_free_symbolic(&Symbolic);
    /* solve system */
    umfpack_di_solve(UMFPACK_A, Ap.ptr, Ai.ptr, Ax.ptr, x.ptr, bx.ptr, Numeric, null, null);
    umfpack_di_free_numeric(&Numeric);
    return x;
}

double[] linearsolve_store_symbolic(SpmCompact!(double) m, out void* Symbolic,
        double[] bx, double[] x = null)
{
    int n = to!int(bx.length);
    void* Numeric;
    if (x is null)
        x = newvec(n, double.init);
    int[] Ap = m.count;
    int[] Ai = m.id;
    double[] Ax = m.data;
    umfpack_di_symbolic(n, n, Ap.ptr, Ai.ptr, Ax.ptr, &Symbolic, null, null);
    umfpack_di_numeric(Ap.ptr, Ai.ptr, Ax.ptr, Symbolic, &Numeric, null, null);
    /* solve system */
    umfpack_di_solve(UMFPACK_A, Ap.ptr, Ai.ptr, Ax.ptr, x.ptr, bx.ptr, Numeric, null, null);
    umfpack_di_free_numeric(&Numeric);
    return x;
}

double[] linearsolve_store_numeric(SpmCompact!(double) m, out void* Numeric,
        double[] bx, double[] x = null)
{
    int n = to!int(bx.length);
    if (x is null)
        x = newvec(n, double.init);
    void* Symbolic;
    int[] Ap = m.count;
    int[] Ai = m.id;
    double[] Ax = m.data;
    umfpack_di_symbolic(n, n, Ap.ptr, Ai.ptr, Ax.ptr, &Symbolic, null, null);
    umfpack_di_numeric(Ap.ptr, Ai.ptr, Ax.ptr, Symbolic, &Numeric, null, null);
    umfpack_di_free_symbolic(&Symbolic);
    /* solve system */
    umfpack_di_solve(UMFPACK_A, Ap.ptr, Ai.ptr, Ax.ptr, x.ptr, bx.ptr, Numeric, null, null);
    return x;
}

double[] linearsolve_with_numeric(SpmCompact!(double) m, void* Numeric,
        double[] bx, double[] x = null)
{
    int n = to!int(bx.length);
    if (x is null)
        x = newvec(n, double.init);
    int[] Ap = m.count;
    int[] Ai = m.id;
    double[] Ax = m.data;
    /* solve system */
    umfpack_di_solve(UMFPACK_A, Ap.ptr, Ai.ptr, Ax.ptr, x.ptr, bx.ptr, Numeric, null, null);
    return x;
}

double[] linearsolve_with_symbolic(SpmCompact!(double) m, void* Symbolic,
        double[] bx, double[] x = null)
{
    int n = to!int(bx.length);
    void* Numeric;
    if (x is null)
        x = newvec(n, double.init);
    int[] Ap = m.count;
    int[] Ai = m.id;
    double[] Ax = m.data;
    umfpack_di_numeric(Ap.ptr, Ai.ptr, Ax.ptr, Symbolic, &Numeric, null, null);
    /* solve system */
    umfpack_di_solve(UMFPACK_A, Ap.ptr, Ai.ptr, Ax.ptr, x.ptr, bx.ptr, Numeric, null, null);
    umfpack_di_free_numeric(&Numeric);
    return x;
}

double[] clinearsolve(CSpmCompact_ri m, double[] b, double[] x = null)
{
    int n = to!int(b.length / 2);
    void* Symbolic, Numeric;
    if (x is null)
        x = newvec(2 * n, double.init);
    int[] Ap = m.count;
    int[] Ai = m.id;
    double[] Ax = m.rdata;
    double[] Az = m.idata;
    double[] xr = x[0 .. n];
    double[] xi = x[n .. 2 * n];
    double[] br = b[0 .. n];
    double[] bi = b[n .. 2 * n];
    umfpack_zi_symbolic(n, n, Ap.ptr, Ai.ptr, Ax.ptr, Az.ptr, &Symbolic, null, null);
    umfpack_zi_numeric(Ap.ptr, Ai.ptr, Ax.ptr, Az.ptr, Symbolic, &Numeric, null, null);
    umfpack_zi_free_symbolic(&Symbolic);
    /* solve system */
    umfpack_zi_solve(UMFPACK_A, Ap.ptr, Ai.ptr, Ax.ptr, Az.ptr, xr.ptr, xi.ptr,
            br.ptr, bi.ptr, Numeric, null, null);
    umfpack_zi_free_numeric(&Numeric);
    return x;
}

double[] linearsolve_store_symbolic(CSpmCompact_ri m, out void* Symbolic,
        double[] b, double[] x = null)
{
    int n = to!int(b.length / 2);
    void* Numeric;
    if (x is null)
        x = newvec(2 * n, double.init);
    int[] Ap = m.count;
    int[] Ai = m.id;
    double[] Ax = m.rdata;
    double[] Az = m.idata;
    double[] xr = x[0 .. n];
    double[] xi = x[n .. 2 * n];
    double[] br = b[0 .. n];
    double[] bi = b[n .. 2 * n];
    umfpack_zi_symbolic(n, n, Ap.ptr, Ai.ptr, Ax.ptr, Az.ptr, &Symbolic, null, null);
    umfpack_zi_numeric(Ap.ptr, Ai.ptr, Ax.ptr, Az.ptr, Symbolic, &Numeric, null, null);
    /* solve system */
    umfpack_zi_solve(UMFPACK_A, Ap.ptr, Ai.ptr, Ax.ptr, Az.ptr, xr.ptr, xi.ptr,
            br.ptr, bi.ptr, Numeric, null, null);
    umfpack_zi_free_numeric(&Numeric);
    return x;
}

double[] linearsolve_store_numeric(CSpmCompact_ri m, out void* Numeric, double[] b, double[] x = null)
{
    int n = to!int(b.length / 2);
    void* Symbolic;
    if (x is null)
        x = newvec(2 * n, double.init);
    int[] Ap = m.count;
    int[] Ai = m.id;
    double[] Ax = m.rdata;
    double[] Az = m.idata;
    double[] xr = x[0 .. n];
    double[] xi = x[n .. 2 * n];
    double[] br = b[0 .. n];
    double[] bi = b[n .. 2 * n];
    umfpack_zi_symbolic(n, n, Ap.ptr, Ai.ptr, Ax.ptr, Az.ptr, &Symbolic, null, null);
    umfpack_zi_numeric(Ap.ptr, Ai.ptr, Ax.ptr, Az.ptr, Symbolic, &Numeric, null, null);
    umfpack_zi_free_symbolic(&Symbolic);
    /* solve system */
    umfpack_zi_solve(UMFPACK_A, Ap.ptr, Ai.ptr, Ax.ptr, Az.ptr, xr.ptr, xi.ptr,
            br.ptr, bi.ptr, Numeric, null, null);
    return x;
}

double[] linearsolve_with_symbolic(CSpmCompact_ri m, void* Symbolic, double[] b, double[] x = null)
{
    int n = to!int(b.length / 2);
    void* Numeric;
    if (x is null)
        x = newvec(2 * n, double.init);
    int[] Ap = m.count;
    int[] Ai = m.id;
    double[] Ax = m.rdata;
    double[] Az = m.idata;
    double[] xr = x[0 .. n];
    double[] xi = x[n .. 2 * n];
    double[] br = b[0 .. n];
    double[] bi = b[n .. 2 * n];
    umfpack_zi_numeric(Ap.ptr, Ai.ptr, Ax.ptr, Az.ptr, Symbolic, &Numeric, null, null);
    /* solve system */
    umfpack_zi_solve(UMFPACK_A, Ap.ptr, Ai.ptr, Ax.ptr, Az.ptr, xr.ptr, xi.ptr,
            br.ptr, bi.ptr, Numeric, null, null);
    umfpack_zi_free_numeric(&Numeric);
    return x;
}

double[] linearsolve_with_numeric(CSpmCompact_ri m, void* Numeric, double[] b, double[] x = null)
{
    int n = to!int(b.length / 2);
    if (x is null)
        x = newvec(2 * n, double.init);
    int[] Ap = m.count;
    int[] Ai = m.id;
    double[] Ax = m.rdata;
    double[] Az = m.idata;
    double[] xr = x[0 .. n];
    double[] xi = x[n .. 2 * n];
    double[] br = b[0 .. n];
    double[] bi = b[n .. 2 * n];
    /* solve system */
    umfpack_zi_solve(UMFPACK_A, Ap.ptr, Ai.ptr, Ax.ptr, Az.ptr, xr.ptr, xi.ptr,
            br.ptr, bi.ptr, Numeric, null, null);
    return x;
}
