/**
 * Basic array and Idx tools.
 *
 * Copyright:   Copyright (C) 2021 by Yinchieh Lai, All Rights Reserved
 * Author:     Yinchieh Lai
 * License:     $(LINK2 http://www.boost.org/LICENSE_1_0.txt, Boost License 1.0)
 */

module addon.tool;
public import std.stdio : write, writeln, writefln, File, readln;
public import std.typetuple: TypeTuple;
public import std.typecons: Tuple, tuple;
public import std.conv: to;
//public import std.traits : CommonType;
import std.process : executeShell;
public import std.format : format;
public import std.parallelism;
public import std.datetime.stopwatch;

alias size_t ST;


public import std.math;
public import std.complex;
alias Complex!float Cfloat;
alias Complex!double Cdouble;
alias Complex!real Creal;
alias complex comp;
alias std.complex.abs cabs;
alias std.complex.arg carg;

import std.variant : Algebraic, VariantN, visit;
alias Algebraic Atype; 

int sys_print(string s)
{
    auto ls = executeShell(s);
    pr(ls.output);
    return ls.status;
}

int sys(string s)
{
    auto ls = executeShell(s);
    if (ls.status != 0)  pr(ls.output);
    return ls.status;
}

/*** templates for array******/

template ArrayBaseType(U)
{
    static if (is(U V == V[]))
    {
        alias ArrayBaseType!(V) ArrayBaseType;
    }
    else
    {
        alias U ArrayBaseType;
    };
}

template NDimArray(U, size_t n)
{
    static if (n == 0)
        alias U NDimArray;
    else
        alias NDimArray!(U, n - 1)[] NDimArray;
}

template ArrayRank(T)
{
    static if (is(T S : S[]))
    {
        const size_t ArrayRank = 1 + ArrayRank!(S);
    }
    else
    {
        const size_t ArrayRank = 0;
    }
}

template ArrayBaseTypeTuple(T...)
{
    static if (T.length == 1)
        alias TypeTuple!(ArrayBaseType!(T[0])) ArrayBaseTypeTuple;
    else
        alias TypeTuple!(ArrayBaseType!(T[0]), ArrayBaseTypeTuple!(T[1 .. $])) ArrayBaseTypeTuple;
}

template NTypeTuple(T, size_t N)
{
    static if (N == 1)
        alias TypeTuple!(T) NTypeTuple;
    else
        alias TypeTuple!(T, NTypeTuple!(T, N - 1)) NTypeTuple;
}

template Iota(size_t a, size_t b)
{
    static if (a < b)
    {
        alias TypeTuple!(a, Iota!(a + 1, b)) Iota;
    }
    else
    {
        alias TypeTuple!() Iota;
    }
}

ref ArrayBaseType!(U) index(U, V...)(U a, V id)
{
    static if (id.length == 1)
        return a[id[0]];
    else
        return index(a[id[0]], id[1 .. $]);
}

template DCType(T, U)
{
    static if (is(T==Cdouble) || is(U==Cdouble))  
        alias Cdouble DCType;
    else
        alias double DCType;
}

/*********function tools for array************/

U[] newvec(U)(size_t dim0, U init = U.init)
{
    U[] sdata = new U[dim0];
    if (init != U.init)
        sdata[] = init;
    return sdata;
}

double[] dvec(size_t dim0, double init = double.init)
{
    return newvec!double(dim0, init);
}

Cdouble[] cvec(size_t dim0, Cdouble init = Cdouble.init)
{
    return newvec!Cdouble(dim0, init);
}

U[] newvec(U)(size_t[] dim, U init = U.init)
{
    return newvec!U(dim[0], init);
}

size_t length_a(U)(U a)
{
    const size_t ndim = ArrayRank!(U);
    static if (ndim == 1)
    {
        return a.length;
    }
    else
    {
        return a.length * length_a(a[0]);
    }
}

ArrayBaseType!(U)* srgptr(U)(U[] a)
{
    const size_t ndim = ArrayRank!(U) + 1;
    static if (ndim == 1)
    {
        return a.ptr;
    }
    else
    {
        return srgptr(a[0]);
    }
}

NDimArray!(ArrayBaseType!(U), ArrayRank!(U) - Ni) subdimarray(size_t Ni, U)(U a)
{
    const size_t ndim = ArrayRank!(U);
    static if (Ni == 0)
    {
        return a;
    }
    else
    {
        return subdimarray!(Ni - 1)(a[0]);
    }
}

ArrayBaseType!(U)[] srg(U)(U[] data)
{
    return (srgptr(data))[0 .. length_a(data)];
}

size_t mod0(V...)(V dim)
{
    static if (V.length == 1)
    {
        return 1;
    }
    else
    {
        return dim[1] * mod0(dim[1 .. $]);
    }
}

auto reshape(W, V...)(W[] sdata0, V dim)
{
    static if (ArrayRank!(W) > 0)
        auto sdata = srg(sdata0);
    else
        auto sdata = sdata0;
    alias ArrayBaseType!(W) U;
    static if (V.length == 1)
        return sdata[0 .. dim[0]];
    else
    {
        size_t d0 = dim[0];
        size_t mod0 = dim.mod0();
        NDimArray!(U, V.length) data = new NDimArray!(U, V.length - 1)[d0];
        size_t itemp = 0;
        foreach (ref v; data)
        {
            v = reshape(sdata[itemp .. itemp + mod0], dim[1 .. $]);
            itemp += mod0;
        }
        return data;
    }
}

auto newarray(U, V...)(V dim)
{
    const size_t n = V.length;
    size_t len = 1;
    foreach (s; dim)
        len *= s;
    return reshape(newvec(len, U.init), dim);
}

auto darray(V...)(V dim)
{
    return newarray!(double)(dim);
}

auto carray(V...)(V dim)
{
    return newarray!(Cdouble)(dim);
}

bool equalp(T)(T[] v1, T[] v2)
{
    if (v1.length != v2.length)
        return false;
    foreach (i; 0 .. v1.length)
    {
        if (v1[i] != v2[i])
            return false;
    }
    return true;
}

size_t[] dims(U)(U u)
{
    const size_t ndim = ArrayRank!(U);
    size_t[] d;
    d.length = ndim;
    foreach (i; Iota!(0, ndim))
    {
        d[i] = subdimarray!(i)(u).length;
    }
    return d;
}

auto dup_a(U)(U[] u)
{
    const size_t ndim = ArrayRank!(U) + 1;
    NTypeTuple!(size_t, ndim) dim;
    foreach (i, ref s; dim)
    {
        s = subdimarray!(i)(u).length;
    }
    return reshape(srg(u), dim);
}

size_t numberOfPoints(T)(T[] x...)
{
    T x0 = x[0];
    T x1 = x[1];
    T dx;
    if (x.length == 3)
        dx = x[2];
    else
        dx = to!T(1);
    size_t np = to!size_t((x1 - x0) / dx) + 1;
    if (x0 + (np - 1) * dx == x1)
        np -= 1;
    return np;
}

T[] arange(T)(T[] x...)
{
    size_t n = numberOfPoints(x);
    T x0 = x[0];
    T x1 = x[1];
    T dx;
    if (x.length == 3)
        dx = x[2];
    else
        dx = to!T(1);
    auto v = newvec!(T)(n);
    foreach (size_t i, ref T s; v)
        s = x0 + i * dx;
    return v;
}

T[] range(T)(T x0, T x1, T dx=cast(T) 1)
{
  return [x0,x1,dx];
}

T[] linspace_a(T)(T x0, T x1, size_t n)
{
    T dx = (x1 - x0) / n;
    auto v = newvec!(T)(n);
    foreach (size_t i, ref T s; v)
        s = x0 +( i+0.5) * dx;
    return v;
}

U arrayeach(alias f, U)(U m)
{
    const size_t ndim = ArrayRank!(U);
    static if (ndim == 1)
        foreach (ref s; m)
            f(s);
    else static if (ndim > 1)
                foreach (v; m)
                    arrayeach!(f)(v);
            else
                assert(false);
    return m;
}

U arrayeach(alias f, U, size_t)(size_t[] id, U m)
{
    const size_t ndim = ArrayRank!(U);
    static if (ndim == 1)
        foreach (i, ref s; m)
            {
            id[$ - ndim] = i;
            f(s);
        }
    else static if (ndim > 1)
        foreach (i, v; m)
            {
            id[$ - ndim] = i;
            arrayeach!(f)(id, v);
        }
    else
        assert(false);
    return m;
}

NDimArray!(V, U.length) table(V, alias f, U...)(U u)
{
    enum n = U.length;
    NTypeTuple!(size_t, n) nv;
    foreach (i, s; u)
        nv[i] = numberOfPoints(s);
    auto m = newarray!(V)(nv);
    size_t[n] id;
    ArrayBaseTypeTuple!(U) var;
    arrayeach!((ref V s) {
        foreach (i, su; u)
        {
            var[i] = su[0] + id[i] * su[2];
        }
        s = f(var);
    })(id, m);
    return m;
}

NDimArray!(double, U.length) dtable(alias f, U...)(U u)
{
    return table!(double, f, U)(u);
}

NDimArray!(Cdouble, U.length) ctable(alias f, U...)(U u)
{
    return table!(Cdouble, f, U)(u);
}

Tuple!(ArrayBaseType!(U), "min", ArrayBaseType!(U), "max") minmax(U)(U[] data)
{
    alias ArrayBaseType!(U) T;
    T max = cast(T)-1.0e20;
    T min = -max;
    foreach (x; srg(data))
    {
        if (x < min)
        {
            min = x;
        }
        if (x > max)
        {
            max = x;
        }
    }
    Tuple!(T, "min", T, "max") result;
    result.min = min;
    result.max = max;
    return result;
}

Tuple!(double, "min", double, "max") absminmax(V)(V[] data)
{
    alias double U;
    alias ArrayBaseType!V T;
    U max = cast(U)-1.0e20;
    U min = -max;
    double temp;
    foreach (x; srg(data))
    {
        temp = abs(x);
        if (temp < min)
        {
            min = temp;
        }
        if (temp > max)
        {
            max = temp;
        }
    }
    Tuple!(U, "min", U, "max") result;
    result.min = min;
    result.max = max;
    return result;
}

double abs2sum(U)(U[] d)
{
    alias ArrayBaseType!U T;
    double temp = 0.;
    foreach (T x; srg(d))
    {
        temp += pow(abs(x), 2);
    }
    return temp;
}

double norm(U)(U[] d)
{
    return sqrt(abs2sum(d));
}

double inner(U, V)(U[] v1, V[] v2)
{
    static if (is(V : Cdouble))
    {
        Cdouble result = comp(0.0, 0.0);
        for (size_t i = 0; i < v1.length; ++i)
            result += v1[i] * conj(v2[i]);
        return result.re;
    }
    else static if (is(U : Cdouble))
    {
        Cdouble result = comp(0.0, 0.0);
        for (size_t i = 0; i < v1.length; ++i)
            result += conj(v1[i]) * v2[i];
        return result.re;
    }
    else
    {
        double result = 0.;
        for (size_t i = 0; i < v1.length; ++i)
            result += v1[i] * v2[i];
        return result;
    }
}

U dot(U)(U[] v1, U[] v2)
{
    U result = cast(U) 0.;
    for (size_t i = 0; i < v1.length; ++i)
        result += v1[i] * v2[i];
    return result;
}

T min(T)(T a, T b)
{
    return (a > b) ? b : a;
}

T max(T)(T a, T b)
{
    return (a > b) ? a : b;
}

NDimArray!(U, ArrayRank!(V[0])) arraymap(U, alias f, W:
        size_t, V...)(W[] id, V v)
{
    enum ndim = ArrayRank!(V[0]);
    NTypeTuple!(size_t, ndim) dim;
    foreach (i, ref s; dim)
    {
        s = subdimarray!(i)(v[0]).length;
    }
    auto m = newarray!(U)(dim);
    ArrayBaseTypeTuple!(V) var;
    NTypeTuple!(size_t, ndim) idtuple;
    arrayeach!((ref U s) {
        foreach (j, ref si; idtuple)
            si = id[j];
        foreach (i, sv; v)
        {
            var[i] = index(sv, idtuple);
        }
        s = f(var);
    })(id, m);
    return m;
}

NDimArray!(U, ArrayRank!(V[0])) arraymap(U, alias f, V...)(V v)
{
    enum ndim = ArrayRank!(V[0]);
    NTypeTuple!(size_t, ndim) dim;
    foreach (i, ref s; dim)
    {
        s = subdimarray!(i)(v[0]).length;
    }
    auto m = newarray!(U)(dim);
    size_t[ndim] id;
    ArrayBaseTypeTuple!(V) var;
    NTypeTuple!(size_t, ndim) idtuple;
    arrayeach!((ref U s) {
        foreach (j, ref si; idtuple)
            si = id[j];
        foreach (i, sv; v)
        {
            var[i] = index(sv, idtuple);
        }
        s = f(var);
    })(id, m);
    return m;
}

W arraymapto(alias f, U:
        size_t, W, V...)(U[] id, W m, V v)
{
    enum ndim = ArrayRank!(W);
    ArrayBaseTypeTuple!(V) var;
    NTypeTuple!(size_t, ndim) idtuple;
    arrayeach!((ref ArrayBaseType!W s) {
        foreach (j, ref si; idtuple)
            si = id[j];
        foreach (i, sv; v)
        {
            var[i] = index(sv, idtuple);
        }
        s = f(var);
    })(id, m);
    return m;
}

W arraymapto(alias f, W, V...)(W m, V v)
{
    enum ndim = ArrayRank!(W);
    size_t[ndim] id;
    ArrayBaseTypeTuple!(V) var;
    NTypeTuple!(size_t, ndim) idtuple;
    arrayeach!((ref ArrayBaseType!W s) {
        foreach (j, ref si; idtuple)
            si = id[j];
        foreach (i, sv; v)
        {
            var[i] = index(sv, idtuple);
        }
        s = f(var);
    })(id, m);
    return m;
}

T[] vec(T, U)(U[] u)
{
 T[] r;
 r.length=u.length;
 foreach(i, s; u) r[i]=to!T(s);
 return r; 
}

T[] vec(T, U...)(U u)
{
 T[] r;
 r.length=U.length;
 foreach(i, s; u) r[i]=to!T(s);
 return r; 
}


/****************generic print function****************************/

void pr1(T)(T t)
{
    static if (is(T == Cdouble) || is(T == Creal))
        write(t, " ");
    else static if (is(T == struct) || is(T == class))
        t.print;
    else static if (isAlgebraic!T)
        t.dispatch!((x){pr(x);});
    else
        write(t, " ");
}

void pr(T...)(T t)
{
    foreach (i, s; t)
    {
        pr1(s);
    }
    writeln("");
}

void pr_obj(T)(T pa)
{
  string[T.tupleof.length] fieldnames; 
  static foreach (i, s; T.tupleof) fieldnames[i]=s.stringof; 
  foreach(i,ref s; pa.tupleof) pr(fieldnames[i]~" = ",s);
}

template Print()
{
  void print()
  {
    string[this.tupleof.length] fieldnames; 
    static foreach (i, s; typeof(this).tupleof) fieldnames[i]=s.stringof;
    pr("struct");
    pr("{");
    foreach(i,ref s; this.tupleof) pr(fieldnames[i]~" = ",s);
    pr("}");
  }
}

/******************file input/output****************/
void arraysave(V)(string fname, V data)
{

    alias ArrayBaseType!(V) U;
    enum ndim = ArrayRank!(V);
    auto fp = File(fname, "w");
    size_t magic;
    static if (is(U == int))
        magic = 0;
    else static if (is(U == long))
        magic = 1;
    else static if (is(U == double))
        magic = 2;
    else static if (is(U == Cdouble))
        magic = 3;
    else static if (is(U == real))
        magic = 4;
    else static if (is(U == Creal) || is(U == creal))
        magic = 5;
    else
        assert(false);
    fp.rawWrite([magic]);
    fp.rawWrite([ndim]);
    fp.rawWrite(dims(data));
    fp.rawWrite(srg(data));
    fp.close();
}

V arrayread(V)(string fname)
{
    alias ArrayBaseType!(V) U;
    enum Ndim = ArrayRank!(V);
    auto fp = File(fname, "r");
    size_t[1] itemp;
    fp.rawRead(itemp);
    size_t magic = itemp[0];
    static if (is(U == int))
        assert(magic == 0);
    else static if (is(U == long))
        assert(magic == 1);
    else static if (is(U == double))
        assert(magic == 2);
    else static if (is(U == Cdouble))
        assert(magic == 3);
    else static if (is(U == real))
        assert(magic == 4);
    else static if (is(U == Creal) || is(U == creal))
        assert(magic == 5);
    else
        assert(false);
    fp.rawRead(itemp);
    size_t ndim = itemp[0];
    assert(ndim == Ndim);
    size_t[] dim;
    dim.length = ndim;
    fp.rawRead(dim);
    NTypeTuple!(size_t, Ndim) dims;
    foreach (i, ref s; dims)
        s = to!(size_t)(dim[i]);
    auto m = newarray!U(dims);
    fp.rawRead(srg(m));
    fp.close();
    return m;
}

/***********template for Idx*******************/

template isIdxType(T)
{
    static if (is(T U == Idx!(U)))
        const bool isIdxType = true;
    else
        const bool isIdxType = false;
}

template IdxBaseType(T)
{
    static if (is(T U == Idx!(U)))
        alias U IdxBaseType;
    else
        alias void IdxBaseType;
}

template IdxBaseTypeTuple(T...)
{
    static if (T.length == 1)
        alias TypeTuple!(IdxBaseType!(T[0])) IdxBaseTypeTuple;
    else
        alias TypeTuple!(IdxBaseType!(T[0]), IdxBaseTypeTuple!(T[1 .. $])) IdxBaseTypeTuple;
}

struct SliceRange
{
    size_t ib;
    size_t ie;
}


/***********Idx struct *******************/

enum _maxdim = 4; 

final struct Idx(T)
{
    size_t ndim;
    size_t[_maxdim] dim;
    size_t[_maxdim] mod;
    size_t offset;
    T[] srg;

    final SliceRange opSlice(size_t n)(size_t i, size_t j)
    {
        auto r = SliceRange(i, j);
        return r;
    }

    final size_t opDollar(size_t i)()
    {
        return dim[i];
    }

    final T opIndex(size_t id0)
    {
        return srg[offset + id0 * mod[0]];
    }

    final T opIndex(size_t id0, size_t id1)
    {
        return srg[offset + id0 * mod[0] + id1 * mod[1]];
    }

    final T opIndex(U : size_t)(U[] args...)
    {
        size_t ip = offset;
        for (size_t i = 0; i < args.length; ++i)
        {
            ip += args[i] * mod[i];
        }
        return srg[ip];
    }

    final Idx!T opIndex(U...)(U u)
    {
        static if (U.length == 0)
            return this.dup();
        else
            return idxslice(this, u);
    }

    final void opIndexAssign(T f)
    {
        if (ndim == 0)
            srg[offset] = f;
        else
            idxeach!((ref x) { x = f; })(this);
    }

    final void opIndexAssign(T[] f)
    {
        assert(this.length() == f.length, "total length not equal");
        size_t id = 0;
        idxeach!((ref x) { x = f[id]; ++id; })(this);
    }

    final void opIndexAssign(Idx!(T) f)
    {
        assert(equaldimp(this, f), "dim not equal");
        size_t[_maxdim] ids;
        size_t[] id = ids[0 .. ndim];
        idxeach_i!((ref x) { x = f[id]; })(id, this);
    }

    final void opIndexAssign(T f, size_t arg0)
    {
        srg[offset + arg0 * mod[0]] = f;
    }

    final void opIndexAssign(T f, size_t arg0, size_t arg1)
    {
        srg[offset + arg0 * mod[0] + arg1 * mod[1]] = f;
    }

    final void opIndexAssign(U : size_t)(T f, U[] args...)
    {
        size_t ip = offset;
        for (size_t i = 0; i < args.length; ++i)
        {
            ip += args[i] * mod[i];
        }
        srg[ip] = f;
    }

    final void opIndexAssign(U...)(T f, U u)
    {
        auto m = idxslice(this, u);
        size_t[_maxdim] ids;
        size_t[] id = ids[0 .. m.ndim];
        idxeach_i!((ref x) { x = f; })(id, m);
    }

    final void opIndexAssign(U...)(T[] f, U u)
    {
        auto m = idxslice(this, u);
        assert(m.length() == f.length, "total length not equal");
        size_t id = 0;
        idxeach!((ref x) { x = f[id]; ++id; })(m);
    }

    final void opIndexAssign(U...)(Idx!(T) f, U u)
    {
        auto m = idxslice(this, u);
        assert(equaldimp(m, f), "dim not equal");
        size_t[_maxdim] ids;
        size_t[] id = ids[0 .. m.ndim];
        idxeach_i!((ref x) { x = f[id]; })(id, m);
    }

    final Idx!T dup()
    {
        auto m0 = newidx!(T)(dim[0 .. ndim]);
        size_t[_maxdim] ids;
        size_t[] id = ids[0 .. ndim];
        idxeach_i!((ref x) { x = this[id]; })(id, m0);
        return m0;
    }

    final int opApply(int delegate(ref size_t, Idx!(T)) dg)
    {
        for (size_t i = 0; i < dim[0]; ++i)
        {
            dg(i, this.bs(i));
        }
        return 0;
    }

    final int opApply(int delegate(ref size_t, ref T) dg)
    {
        assert(ndim == 1);
        size_t inc = mod[0];
        T* p = srg.ptr + offset;
        for (size_t i = 0; i < dim[0]; ++i)
        {
            dg(i, *p);
            p += inc;
        }
        return 0;
    }

    final int opApply(int delegate(Idx!(T)) dg)
    {
        for (size_t i = 0; i < dim[0]; ++i)
        {
            dg(this.bs(i));
        }
        return 0;
    }

    final int opApply(int delegate(ref T) dg)
    {
        assert(ndim == 1);
        T* p = srg.ptr + offset;
        size_t inc = mod[0];
        for (size_t i = 0; i < dim[0]; ++i)
        {
            dg(*p);
            p += inc;
        }
        return 0;
    }

    final int opApplyReverse(int delegate(ref size_t, Idx!(T)) dg)
    {
        for (size_t i = 0; i < dim[0]; ++i)
        {
            dg(i, this.es(i));
        }
        return 0;
    }

    final int opApplyReverse(int delegate(Idx!(T)) dg)
    {
        for (size_t i = 0; i < dim[0]; i++)
        {
            dg(this.es(i));
        }
        return 0;
    }

    /****methods *******************/

    final size_t length()
    {
        assert(ndim > 0, "wrong ndim!");
        size_t result = 1;
        for (size_t i = 0; i < ndim; ++i)
        {
            result *= dim[i];
        }
        return result;
    }

    final size_t[] dims()
    {
        assert(ndim > 0, "wrong ndim!");
        return dim[0 .. ndim];
    }

    final auto bs(size_t i0)
    {
        Idx!(T) v;
        v.srg = srg;
        v.ndim = ndim - 1;
        v.offset = offset + i0 * mod[0];
        for (size_t i = 1; i < ndim; ++i)
        {
            v.dim[i - 1] = dim[i];
            v.mod[i - 1] = mod[i];
        }
        return v;
    }

    final auto es(size_t i0)
    {
        Idx!(T) v;
        v.srg = srg;
        v.ndim = ndim - 1;
        v.offset = offset + i0 * mod[ndim - 1];
        for (size_t i = 0; i < ndim - 1; ++i)
        {
            v.dim[i] = dim[i];
            v.mod[i] = mod[i];
        }
        return v;
    }

    final NDimArray!(T, N) srg_a(size_t N)()
    {
        assert(this.contiguousp, "idx not contiguous");
        NTypeTuple!(size_t, N) var;
        foreach (i, ref s; var)
        {
            s = this.dim[i];
        }
        return reshape(this.srg, var);
    }

    final Idx!(T) sampling(size_t sdim = 0)(size_t ns)
    {
        return resize(this, sdim, [0, dim[sdim], ns]);
    }

    final void print()
    {
        idxprint(this);
    }

    final Idx!(T) opAssign(Idx!T m0)
    {
        ndim = m0.ndim;
        dim[0 .. ndim] = m0.dim[0 .. ndim];
        mod[0 .. ndim] = m0.mod[0 .. ndim];
        offset = m0.offset;
        srg = m0.srg;
        return this;
    }

}

/***********function tools for Idx*************************/

alias Idx!(double) DIdx;
alias Idx!(Cdouble) CIdx;
alias Idx!(int) IIdx;

Idx!(T) srgtoidx(T)(T[] srg, size_t[] d)
{
    Idx!(T) m0;
    m0.ndim = d.length;
    size_t ilength = 1;
    for (size_t i = m0.ndim - 1; i < m0.ndim; --i)
    {
        m0.dim[i] = d[i];
        m0.mod[i] = ilength;
        ilength *= d[i];
    }
    m0.offset = 0;
    m0.srg = srg;
    return m0;
}

Idx!(T) reshape(T)(Idx!(T) init, size_t[] d)
{
    assert(contiguousp(init), "idx not contiguous");
    return srgtoidx(init.srg, d);
}

Idx!(ArrayBaseType!(U)) idx(U)(U[] data)
{
    return srgtoidx(srg(data), dims(data));
}

Idx!(ArrayBaseType!(U)) idx(U)(U[] data, size_t[] d)
{
    return srgtoidx(srg(data), d);
}

Idx!(T) newidx(T, U:
        size_t)(U[] dim, T init = T.init)
{
    size_t len = 1;
    foreach (s; dim)
    {
        len *= s;
    }
    return srgtoidx(newvec!T(len, init), dim);
}

auto newidx_(U, V...)(V dim)
{
    const size_t n = V.length;
    size_t[n] dims;
    size_t len = 1;
    foreach (i, s; dim)
    {
        len *= s;
        dims[i] = s;
    }
    return srgtoidx(newvec!U(len), dims);
}

Idx!(double) didx(size_t[] d, double init = double.init)
{
    return newidx!(double)(d, init);
}

Idx!(Cdouble) cidx(size_t[] d, Cdouble init = Cdouble.init)
{
    return newidx!(Cdouble)(d, init);
}

DIdx didx(V...)(V dim)
{
    return newidx_!(double)(dim);
}

Idx!U emptydup(U, T)(Idx!T m)
{
    return newidx!(U)(m.dims());
}

CIdx cidx(V...)(V dim)
{
    return newidx_!(Cdouble)(dim);
}

template IdxNull(T)
{
    Idx!(T) IdxNull;
}

alias IdxNull!(double) DIdxNull;
alias IdxNull!(Cdouble) CIdxNull;
alias IdxNull!(int) IIdxNull;

bool equaldimp(U, V)(Idx!(U) m1, Idx!(V) m2)
{
    return (m1.ndim == m2.ndim) && (m1.dim[0 .. m1.ndim].equalp(m2.dim[0 .. m1.ndim]));
}

bool contiguousp(U)(Idx!(U) m)
{
    size_t len = 1;
    for (size_t i = m.ndim - 1; i < m.ndim; --i)
    {
        if (m.mod[i] != len)
        {
            return false;
        }
        len *= m.dim[i];
    }
    return (len == m.srg.length);
}

Idx!(T) slice(T)(Idx!(T) m0, size_t dimension, size_t id)
{
    auto m = m0;
    size_t temp = m.mod[dimension];
    for (size_t i = dimension + 1; i < m.ndim; ++i)
    {
        m.dim[i - 1] = m.dim[i];
        m.mod[i - 1] = m.mod[i];
    }
    m.offset += id * temp;
    m.ndim -= 1;
    return m;
}

Idx!(T) resize(T)(Idx!(T) m0, size_t n, size_t[] arg)
{
    assert((n < m0.ndim) && (arg[0] < m0.dim[n]) && (arg[1] <= m0.dim[n]));
    auto m = m0;
    if (arg.length == 3)
    {
        if ((arg[1] - arg[0]) % arg[2] == 0)
            m.dim[n] = (arg[1] - arg[0]) / arg[2];
        else
            m.dim[n] = (arg[1] - arg[0]) / arg[2] + 1;
        m.mod[n] = arg[2] * m.mod[n];
        m.offset += arg[0] * m.mod[n];
    }
    else if (arg.length == 2)
    {
        m.dim[n] = (arg[1] - arg[0]);
        m.offset += arg[0] * m.mod[n];
    }
    else
        assert(false, "number of index mismatch");
    return m;
}

Idx!(T) idxslice(T, U...)(Idx!T m, U u)
{
    static if (U.length == 0)
        return m;
    else static if (is(U[$ - 1] == SliceRange))
        return idxslice(resize(m, u.length - 1, [u[$ - 1].ib, u[$ - 1].ie]), u[0 .. $ - 1]);
    else
        return idxslice(slice(m, u.length - 1, u[$ - 1]), u[0 .. $ - 1]);
}

V[0] idxeach_i(alias f, U:
        size_t, V...)(U[] id, V m)
{
    size_t ndim = m[0].ndim;
    if (ndim == 1)
        foreach (k, ref IdxBaseType!(V[0]) s0; m[0])
        {
            id[$ - ndim] = k;
            IdxBaseTypeTuple!V var_;
            auto var = var_[1 .. $];
            foreach (i, ref s; var)
                s = m[i + 1][k];
            f(s0, var);
        }
    else if (ndim > 1)
        foreach (k, V[0] v0; m[0])
        {
            id[$ - ndim] = k;
            V var_;
            auto var = var_[1 .. $];
            foreach (i, ref s; var)
                s = m[i + 1].bs(k);
            idxeach_i!(f)(id, v0, var);
        }
    return m[0];
}

V[0] idxeach(alias f, V...)(V m)
{
    size_t ndim = m[0].ndim;
    if (ndim == 1)
        foreach (k, ref IdxBaseType!(V[0]) s0; m[0])
        {
            IdxBaseTypeTuple!V var_;
            auto var = var_[1 .. $];
            foreach (i, ref s; var)
                s = m[i + 1][k];
            f(s0, var);
        }
    else if (ndim > 1)
        foreach (k, V[0] v0; m[0])
        {
            V var_;
            auto var = var_[1 .. $];
            foreach (i, ref s; var)
                s = m[i + 1].bs(k);
            idxeach!(f)(v0, var);
        }
    return m[0];
}

Idx!(V) idxtable(V, alias ft, U...)(U u)
{
    enum n = U.length;
    size_t[n] nv;
    foreach (i, su; u)
        nv[i] = numberOfPoints(su);
    auto m = newidx!(V)(nv);
    size_t[n] id;
    ArrayBaseTypeTuple!(U) var;
    idxeach_i!((ref s) {
        foreach (i, su; u)
        {
            var[i] = su[0] + id[i] * su[2];
        }
        s = ft(var);
    })(id, m);
    return m;
}

Idx!(double) didxtable(alias f, U...)(U u)
{
    return idxtable!(double, f, U)(u);
}

Idx!(Cdouble) cidxtable(alias f, U...)(U u)
{
    return idxtable!(Cdouble, f, U)(u);
}

Tuple!(T, "min", T, "max") minmax(T)(Idx!(T) data)
{
    T max = cast(T)-1.0e20;
    T min = -max;
    T temp;
    idxeach!((T x) {
        temp = x;
        if (temp < min)
        {
            min = temp;
        }
        if (temp > max)
        {
            max = temp;
        }
    })(data);
    Tuple!(T, "min", T, "max") result;
    result.min = min;
    result.max = max;
    return result;
}

Tuple!(double, "min", double, "max") absminmax(T)(Idx!(T) data)
{
    double max = -1.0e20;
    double min = -max;
    double temp;
    idxeach!((T x) {
        temp = abs(x);
        if (temp < min)
        {
            min = temp;
        }
        if (temp > max)
        {
            max = temp;
        }
    })(data);
    Tuple!(double, "min", double, "max") result;
    result.min = min;
    result.max = max;
    return result;
}

double abs2sum(T)(Idx!(T) d)
{
    double temp = 0.;
    idxeach!((T x) { temp += pow(abs(x), 2); })(d);
    return temp;
}

Idx!T idxrange(T)(T[] x...)
{
    return arange(x).idx;
}

double norm(T)(Idx!(T) d)
{
    return sqrt(abs2sum(d));
}

Idx!(U) join(size_t n, U)(Idx!(U) m1, Idx!(U) m2)
{
    assert(m1.ndim == m2.ndim, "dim error");
    assert(m1.dim[0 .. n].equalp(m2.dim[0 .. n])
            && m1.dim[n + 1 .. m1.ndim].equalp(m2.dim[n + 1 .. m2.ndim]), "dim error");
    size_t[] dim = m1.dims.dup;
    dim[n] += m2.dim[n];
    auto m = newidx!(U)(dim);
    m.resize(n, [0, m1.dim[n]])[] = m1;
    m.resize(n, [m1.dim[n], m.dim[n]])[] = m2;
    return m;
}

void idxprint(T)(Idx!(T) m)
{
    if (m.ndim == 1)
    {
        write("{");
        foreach (T s; m)
            write(s, " ");
        writeln("}");
    }
    else if (m.ndim > 1)
    {
        write("{");
        foreach (Idx!(T) v; m)
        {
            idxprint(v);
        }
        writeln("}");
    }
    else if (m.ndim == 0)
        writeln("{", m.srg[m.offset], "} ");
    else
        assert(false, "ndim error");
}


DIdx re(CIdx m)
{
    assert(m.contiguousp(), "idx not contiguous.");
    size_t[_maxdim] newdim;
    newdim[0 .. m.ndim][] = m.dims;
    newdim[m.ndim] = 2;
    return (cast(double*) m.srg.ptr)[0 .. 2 * m.length].idx(newdim[0 .. m.ndim + 1]).es(0);
}

DIdx im(CIdx m)
{
    assert(m.contiguousp(), "idx not contiguous.");
    size_t[_maxdim] newdim;
    newdim[0 .. m.ndim][] = m.dims;
    newdim[m.ndim] = 2;
    return (cast(double*) m.srg.ptr)[0 .. 2 * m.length].idx(newdim[0 .. m.ndim + 1]).es(1);
}

version (LDC)
{

}
else
{
    Idx!U idxmap_i(U, alias f, W, V...)(W[] id, V m)
    {
        return idxeach_i!((ref U s, IdxBaseTypeTuple!(V) var) { s = f(var); })(id,
                m[0].emptydup!U, m);
    }

    Idx!U idxmapto_i(alias f, W:
            size_t, U, V...)(W[] id, Idx!U m0, V m)
    {
        return idxeach_i!((ref U s, IdxBaseTypeTuple!(V) var) { s = f(var); })(id, m0, m);
    }
}

Idx!U idxmap(U, alias f, V...)(V m)
{
    return idxeach!((ref U s, IdxBaseTypeTuple!(V) var) { s = f(var); })(m[0].emptydup!U, m);
}

Idx!U idxmapto(alias f, U, V...)(Idx!U m0, V m)
{
    return idxeach!((ref U s, IdxBaseTypeTuple!(V) var) { s = f(var); })(m0, m);
}

CIdx to_cidx(DIdx m0)
{
  auto m=newidx!(Cdouble)(m0.dims);
  size_t[] ip;
  ip.length=m0.ndim;
  idxeach_i!((ref s){ s=comp(m0[ip]);})(ip,m);
  return m;  
}

DIdx cabs(CIdx m)
{
  return idxmap!(double, (z)=>cabs(z))(m);
}

DIdx cabs2(CIdx m)
{
  return idxmap!(double, (z)=>cabs(z)^^2)(m);
}

DIdx carg(CIdx m)
{
  return idxmap!(double, (z)=>carg(z))(m);
}


Idx!(T) transpose(T)(Idx!(T) m0)
{
    if (m0.ndim == 2)
    {
    auto m = m0;
    m.dim[1] = m0.dim[0];
    m.dim[0] = m0.dim[1];
    m.mod[1] = m0.mod[0];
    m.mod[0] = m0.mod[1];
    return m;
  }
  else if (m0.ndim ==1) return m0;
  else assert(false, "not implemented");
}


/******************************/
template declare_incr(u...)
{
    static if (u.length == 2)
        const string declare_incr = format(r"size_t %s_incr=%s.mod[0];
    ", u[1], u[1]);
    else
        const string declare_incr = format(r"size_t %s_incr=%s.mod[0];
  ", u[1], u[1]) ~ declare_incr!(u[2 .. $]);
}

template declare_pointer(u...)
{
    static if (u.length == 2)
        const string declare_pointer = format(r"auto %s = &( %s.srg[ %s.offset]);
    ",
                u[0], u[1], u[1]);
    else
        const string declare_pointer = format(r"auto %s = &( %s.srg[ %s.offset]);
  ",
                u[0], u[1], u[1]) ~ declare_pointer!(u[2 .. $]);
}

template update_pointer(u...)
{
    static if (u.length == 2)
        const string update_pointer = format(r"%s += %s_incr;
    ", u[0], u[1]);
    else
        const string update_pointer = format(r"%s += %s_incr;
  ", u[0], u[1]) ~ update_pointer!(u[2 .. $]);
}

template forloop_begin(string u)
{
    const string forloop_begin = format(r"for(size_t i=0;i < %s.dim[0];++i){ 
   ", u);
}

template idx1each(U...)
{
    const string idx1each = r"{" ~ assert_idx1!(U[0 .. $ - 1]) ~ declare_incr!(
            U[0 .. $ - 1]) ~ declare_pointer!(U[0 .. $ - 1]) ~ forloop_begin!(
            U[1]) ~ U[$ - 1] ~ update_pointer!(U[0 .. $ - 1]) ~ r"}" ~ r"}";
}

template assert_idx1(u...)
{
    static if (u.length == 2)
        const string assert_idx1 = format("assert(%s.ndim==1,  \"ndim error\");", u[1]);
    else
        const string assert_idx1 = format("assert(%s.ndim==1,  \"ndim error\");", u[1])
            ~ assert_idx1!(u[2 .. $]);
}

DCType!(T, U) idx1dot1(T, U)(Idx!(T) v1, Idx!(U) v2)
{
    assert(equaldimp(v1, v2), "unmatched dim");
    DCType!(T, U) result = cast(DCType!(T, U)) 0.;
    mixin(idx1each!("p1", "v1", "p2", "v2", "result += (*p1)*(*p2);"));
    return result;
}

DCType!(T, U) idx1cdot1(T, U)(Idx!(T) v1, Idx!(U) v2)
{
    assert(equaldimp(v1, v2), "unmatched dim");
    DCType!(T, U) result = cast(DCType!(T, U)) 0.;
    mixin(idx1each!("p1", "v1", "p2", "v2", "result += conj(*p1)*(*p2);"));
    return result;
}

DCType!(T, U) idx2dot1(T, U, V = DCType!(T, U))(Idx!(T) ma,
        Idx!(U) vb, Idx!(V) vc = IdxNull!(V))
{
    assert(ma.ndim == 2 && vb.ndim == 1);
    if (vc.ndim == 0)
        vc = newidx_!V([ma.dim[0]]);
    foreach (size_t i, Idx!(T) va; ma)
    {
        vc[i] = idx1dot1(va, vb);
    }
    return vc;
}

DCType!(T, U) idx1dot2(T, U, V =DCType!(T, U))(Idx!(T) va,
        Idx!(U) mb, Idx!(V) vc = IdxNull!(V))
{
    assert(va.ndim == 1 && mb.ndim == 2);
    if (vc.ndim == 0)
        vc = newidx_!V([mb.dim[1]]);
    foreach_reverse (size_t i, Idx!(T) vb; mb)
    {
        vc[i] = idx1dot1(va, vb);
    }
    return vc;
}

DCType!(T, U) idx2dot2(T, U, V = DCType!(T, U))(Idx!(T) m1,
        Idx!(U) m2, Idx!(V) m = IdxNull!(V))
{
    assert(m1.ndim == 2 && m2.ndim == 2);
    if (m.ndim == 0)
        m = newidx_!V([m1.dim[0], m2.dim[1]]);
    foreach_reverse (size_t i, Idx!(V) v; m)
        idx2dot1(m1, m2.es(i), v);
    return m;
}

double inner(U, V)(Idx!(U) v1, Idx!(V) v2)
{
    assert(equaldimp(v1, v2), "unmatched dim");
    static if (is(U : Cdouble) || is(V : Cdouble))
    {
        Cdouble result = comp(0.0, 0.0);
        mixin(idx1each!("p1", "v1", "p2", "v2", "result += conj(*p1)*(*p2);"));
        return result.re;
    }
    else
    {
        double result = 0.;
        mixin(idx1each!("p1", "v1", "p2", "v2", "result += (*p1)*(*p2);"));
        return result;
    }
}

/******************Idx file input/output****************/

void idxsave(U)(string fname, Idx!(U) data)
{
    assert(data.contiguousp(), "idx not contiguous.");
    auto fp = File(fname, "w");
    int magic;
    static if (is(U == int))
        magic = 0;
    else static if (is(U == long))
        magic = 1;
    else static if (is(U == double))
        magic = 2;
    else static if (is(U == Cdouble))
        magic = 3;
    else static if (is(U == real))
        magic = 4;
    else static if (is(U == Creal) || is(U == creal))
        magic = 5;
    else
        assert(false);
    fp.rawWrite([magic]);
    fp.rawWrite([to!int(data.ndim)]);
    fp.rawWrite(to!(int[])(data.dim[0 .. data.ndim]));
    fp.rawWrite(data.srg);
    fp.close();
}

Idx!(U) idxread(U)(string fname)
{
    auto fp = File(fname, "r");
    int[1] itemp;
    fp.rawRead(itemp);
    int magic = itemp[0];
    static if (is(U == int))
        assert(magic == 0);
    else static if (is(U == long))
        assert(magic == 1);
    else static if (is(U == double))
        assert(magic == 2);
    else static if (is(U == Cdouble))
        assert(magic == 3);
    else static if (is(U == real))
        assert(magic == 4);
    else static if (is(U == Creal) || is(U == creal))
        assert(magic == 5);
    else
        assert(false);
    fp.rawRead(itemp);
    size_t ndim = to!(size_t)(itemp[0]);
    int[] dim;
    dim.length = ndim;
    fp.rawRead(dim);
    auto m = newidx!U(to!(size_t[])(dim));
    fp.rawRead(m.srg);
    fp.close();
    return m;
}

/**********grid points in polygon********************/
pragma(inline) long unitstep(double x)
{
  return (x<0)?0:1;
}

void rotateleft(T)(T[] u, T[] r)
{
  size_t np=u.length;
  for(long i=0;i<np-1;++i)
  {
    r[i]=u[i+1]; 
  }
  r[np-1]=u[0]; 
}

void setPoly(T)(double[][] poly, Idx!T m, T val, double[] x, double[] y)
{
  long np=poly.length;
  double[] Xi,Yi,Xip1,Yip1;
  long[] u,v,w;
  Xi.length=Yi.length=Xip1.length=Yip1.length=u.length=v.length=w.length=np;

  bool inPolyQ(double x, double y)
  {
    for(long i=0;i<np;++i)
    {
      Xi[i]=poly[i][0]; 
      Yi[i]=poly[i][1]; 
    }
    for(long i=0;i<np-1;++i)
    {
      Xip1[i]=Xi[i+1]; 
      Yip1[i]=Yi[i+1]; 
    }
    Xip1[np-1]=Xi[0]; 
    Yip1[np-1]=Yi[0]; 
    for(long i=0;i<np;++i)
    {
      u[i]=unitstep(y-Yi[i]);  
    }
    rotateleft(u,v);
    for(long i=0;i<np;++i)
    {
      w[i]=unitstep(-((Xip1[i] - Xi[i])*(y - Yi[i]) - (x - Xi[i])*(Yip1[i] - Yi[i])));  
    }
    long sum=0;
    for(long i=0;i<np;++i)
    {
      sum += (u[i]*(1 - v[i])*(1 - w[i]) - (1 - u[i])*v[i]*w[i]);  
    }
    return (sum != 0)?true:false;
  }

  ST[2] ip;
  idxeach_i!(
    (ref s){
      if(inPolyQ(x[ip[0]],y[ip[1]])) s=val;
    }) (ip,m);
  Xi.destroy();
  Yi.destroy();
  Xip1.destroy();
  Yip1.destroy();
  u.destroy();
  v.destroy();
  w.destroy();
}

/*********templaye for Algebraic type ***************************/

template isAlgebraic(Type)
{
    static if (is(Type == VariantN!T, T...))
        enum isAlgebraic = T.length >= 2;
    else
        enum isAlgebraic = false;
}

template dispatch(alias Fn)
{
    void dispatch(A)(A arg)
    if (A.AllowedTypes.length > 0)
    {
        static foreach (T; A.AllowedTypes)
        {
            if (typeid(T) == arg.type) {
                Fn(arg.get!T);
            }
        }
    }
}

/+ example:
  alias Atype!(Cdouble, int, double) Mytype;
  auto v=vec!Mytype(1, comp(2.), 3.);
  pr(v);
+/
