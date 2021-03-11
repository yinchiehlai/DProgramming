/**
 * Basic vector, matrix, idx tools.
 *
 * Copyright:   Copyright (C) 2021 by Yinchieh Lai
 * Author:     Yinchieh Lai
 * License:     $(LINK2 http://www.boost.org/LICENSE_1_0.txt, Boost License 1.0)
 */

module addon.spm;
import addon.tool;

struct Spme(T)
{
    size_t id;
    T data;
    Spme!(T)* cdr;

    final int opApply(int delegate(Spme!(T)) dg)
    {
        Spme!(T)* p = &this;
        while (p !is null)
        {
            dg(*p);
            p = p.cdr;
        }
        return 0;
    }

    final size_t length()
    {
        size_t n = 0;
        foreach (Spme!(T) s; this)
        {
            ++n;
        }
        return n;
    }
}

/*************************************************/
Spme!(T)* spmcons(T)(size_t id, T data, Spme!(T)* cdr, Spme!(T)* p)
{
    p.id = id;
    p.data = data;
    p.cdr = cdr;
    return p;
}

void spm_scan(alias f, T)(Spm!T m)
{
    for (size_t j = 0; j < m.dim[1]; ++j)
    {
        if (m.col[j]!is null)
        {
            foreach (s; *(m.col[j]))
                f(s.id, j, s.data);
        }
    }
}

struct Spm(T)
{
    size_t[2] dim;
    T def;
    Spme!(T)*[] col;
    size_t len;
    Spme!(T)[] buf;
    size_t ib;

    final T opIndex(size_t i, size_t j)
    {
        assert(j < dim[1], "Spm dim1 bound error");
        auto p = col[j];
        size_t ip;
        for (;;)
        {
            if (p is null)
            {
                return def;
            }
            else
            {
                ip = p.id;
                if (i < ip)
                    return def;
                else if (i == ip)
                    return p.data;
                else
                    p = p.cdr;
            }
        }
    }

    final T opIndexAssign(T data, size_t i, size_t j)
    {
        assert(j < dim[1], "Spm dim1 bound error");
        auto p = col[j];
        size_t ip;
        if (p is null)
        {
            col[j] = spmcons(i, data, null, &buf[ib++]);
            ++len;
        }
        else
        {
            ip = p.id;
            if (i < ip)
            {
                p.cdr = spmcons(p.id, p.data, p.cdr, &buf[ib++]);
                p.id = i;
                p.data = data;
                ++len;
            }
            else if (i == ip)
            {
                p.data = data;
            }
            else
            {
                while (p.cdr !is null)
                {
                    ip = p.cdr.id;
                    if (i < ip)
                    {
                        p.cdr = spmcons(i, data, p.cdr, &buf[ib++]);
                        ++len;
                        return data;
                    }
                    else if (i == ip)
                    {
                        p.cdr.data = data;
                        return data;
                    }
                    else
                        p = p.cdr;
                }
                p.cdr = spmcons(i, data, null, &buf[ib++]);
                ++len;
            }
        }
        return data;
    }

    final int opApply(int delegate(Spme!(T)) dg)
    {
        foreach (Spme!(T)* v; col)
        {
            if (v !is null)
                foreach (Spme!(T) s; *v)
                {
                    dg(s);
                }
        }
        return 0;
    }

    final int opApply(int delegate(ref size_t, Spme!(T)) dg)
    {
        foreach (j, Spme!(T)* v; col)
        {
            if (v !is null)
                foreach (Spme!(T) s; *v)
                {
                    dg(j, s);
                }
        }
        return 0;
    }

    final void opIndexOpAssign(string s0)(T data, size_t i, size_t j)
    {
        static if (s0 == "+")
        {
            assert(j < dim[1], "Spm dim1 bound error");
            auto p = col[j];
            size_t ip;
            if (p is null)
            {
                col[j] = spmcons(i, data, null, &buf[ib++]);
                ++len;
            }
            else
            {
                ip = p.id;
                if (i < ip)
                {
                    p.cdr = spmcons(p.id, p.data, p.cdr, &buf[ib++]);
                    p.id = i;
                    p.data = data;
                    ++len;
                }
                else if (i == ip)
                {
                    p.data += data;
                }
                else
                {
                    while (p.cdr !is null)
                    {
                        ip = p.cdr.id;
                        if (i < ip)
                        {
                            p.cdr = spmcons(i, data, p.cdr, &buf[ib++]);
                            ++len;
                            return;
                        }
                        else if (i == ip)
                        {
                            p.cdr.data += data;
                            return;
                        }
                        else
                            p = p.cdr;
                    }
                    p.cdr = spmcons(i, data, null, &buf[ib++]);
                    ++len;
                }
            }
        }
        else static if (s0 == "-")
        {
            assert(j < dim[1], "Spm dim1 bound error");
            auto p = col[j];
            size_t ip;
            if (p is null)
            {
                col[j] = spmcons(i, -data, null, &buf[ib++]);
                ++len;
            }
            else
            {
                ip = p.id;
                if (i < ip)
                {
                    p.cdr = spmcons(p.id, p.data, p.cdr, &buf[ib++]);
                    p.id = i;
                    p.data = -data;
                    ++len;
                }
                else if (i == ip)
                {
                    p.data -= data;
                }
                else
                {
                    while (p.cdr !is null)
                    {
                        ip = p.cdr.id;
                        if (i < ip)
                        {
                            p.cdr = spmcons(i, -data, p.cdr, &buf[ib++]);
                            ++len;
                            return;
                        }
                        else if (i == ip)
                        {
                            p.cdr.data -= data;
                            return;
                        }
                        else
                            p = p.cdr;
                    }
                    p.cdr = spmcons(i, -data, null, &buf[ib++]);
                    ++len;
                }
            }
        }
    }

    final CSpmCompact_ri compact_ri()
    {
        CSpmCompact_ri r;
        r.count = newvec!(int)(dim[1] + 1);
        r.count[0] = 0;
        for (size_t j = 1; j <= dim[1]; ++j)
        {
            if (col[j - 1]!is null)
                r.count[j] = cast(int) col[j - 1].length + r.count[j - 1];
            else
                r.count[j] = r.count[j - 1];
        }
        size_t np = r.count[dim[0]];
        r.id = newvec!(int)(np);
        r.rdata = newvec!(double)(np);
        r.idata = newvec!(double)(np);
        r.diagonal = newvec!(size_t)(np);
        size_t k = 0;
        foreach (j, Spme!(T) s; this)
        {
            version (X86_64)
                r.id[k] = cast(int) s.id;
            else
                r.id[k] = s.id;
            r.rdata[k] = s.data.re;
            r.idata[k] = s.data.im;
            if (s.id == j)
                r.diagonal[j] = k;
            ++k;
        }
        return r;
    }

    final SpmCompact!T compact()
    {
        SpmCompact!T r;
        r.count = newvec!(int)(dim[1] + 1);
        r.count[0] = 0;
        for (size_t j = 1; j <= dim[1]; ++j)
        {
            if (col[j - 1]!is null)
                r.count[j] = cast(int) col[j - 1].length + r.count[j - 1];
            else
                r.count[j] = r.count[j - 1];
        }
        size_t np = r.count[dim[0]];
        r.id = newvec!(int)(np);
        r.data = newvec!(T)(np);
        r.diagonal = newvec!(size_t)(np);
        size_t k = 0;
        foreach (j, Spme!(T) s; this)
        {
            version (X86_64)
                r.id[k] = cast(int) s.id;
            else
                r.id[k] = s.id;
            r.data[k] = s.data;
            if (s.id == j)
                r.diagonal[j] = k;
            ++k;
        }
        return r;
    }

    final Spm!T dup()
    {
        Spm!T m = newspm!(T)(dim, 0);
        m.buf = newvec!(Spme!T)(buf.length);
        spm_scan!((i, j, data) { m[i, j] = data; })(this);
        return m;
    }

    ~this()
    {
        col.destroy();
        buf.destroy();
    }

    final void print()
    {
        spm_scan!((i, j, data) { pr(i, ",", j, ":", data); })(this);
    }
}

Spm!(T) newspm(T, U:
        size_t)(U[] dim, size_t bufsize, T init = cast(T) 0.)
{
    Spm!(T) m;
    m.dim[0] = dim[0];
    m.dim[1] = dim[1];
    m.col = newvec!(Spme!T*)(dim[1]);
    foreach (ref s; m.col)
        s = null;
    m.len = 0;
    m.def = init;
    if (bufsize != 0)
        m.buf = newvec!(Spme!T)(bufsize);
    m.ib = 0;
    return m;
}

alias Spme!double DSpme;
alias Spme!Cdouble CSpme;
alias Spm!double DSpm;
alias Spm!Cdouble CSpm;

DSpm dspm(U : size_t)(U[] dim, size_t bufsize)
{
    return newspm!(double)(dim, bufsize);
}

CSpm cspm(U : size_t)(U[] dim, size_t bufsize)
{
    return newspm!(Cdouble)(dim, bufsize);
}

struct CSpmCompact_ri
{
    int[] count;
    int[] id;
    double[] rdata;
    double[] idata;
    size_t[] diagonal;

    ~this()
    {
        count.destroy();
        id.destroy();
        rdata.destroy();
        idata.destroy();
        diagonal.destroy();
    }
}

struct SpmCompact(T)
{
    int[] count;
    int[] id;
    T[] data;
    size_t[] diagonal;

    ~this()
    {
        count.destroy();
        id.destroy();
        data.destroy();
        diagonal.destroy();
    }
}

alias SpmCompact!double DSpmCompact;
alias SpmCompact!Cdouble CSpmCompact;

/////////////////////////////////////////
