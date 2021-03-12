/**
 * Basic plotting tools based on the DISLIN library (https://www.mps.mpg.de/dislin).
 *
 * Copyright:   Copyright (C) 2021 by Yinchieh Lai
 * Author:     Yinchieh Lai
 * License:     $(LINK2 http://www.boost.org/LICENSE_1_0.txt, Boost License 1.0)
 */

module addon.plot;
import addon.dislin;
import addon.tool;
import std.file;

pragma(lib, "dislin_d");
version (Windows)
{
    pragma(lib, "gdi32");
    pragma(lib, "user32");
}
else
{
    pragma(lib, "Xm");
}

void listplot(U...)(U u)
{
    movieplot_begin();
    static if (U.length == 1)
    {
        auto x = arange(0., to!double(u[0].length), 1.);
        listplot_m(u[0], x);
        x.destroy();
    }
    else static if (U.length == 2)
    {
        listplot_m(u[0], u[1]);
    }
    plot_end();
}

void listplot_file(U...)(string filenamenoext, U u)
{
    string filename=filenamenoext ~ ".png";
    if (filename.exists) remove(filename);
    paperplot_begin(filename);
    static if (U.length == 1)
    {
        auto x = arange(0., to!double(u[0].length), 1.);
        listplot_m(u[0], x);
        x.destroy();
    }
    else static if (U.length == 2)
    {
        listplot_m(u[0], u[1]);
    }
    plot_end();
}

void listplot3d(double[][] data, double[] xray = null, double[] yray = null)
{
    auto dim = data.dims();
    if (xray is null)
        xray = arange(0., to!(double)(dim[0]), 1.);
    if (yray is null)
        yray = arange(0., to!(double)(dim[1]), 1.);
    auto temp = minmax(data);
    qplsur(data.srg.ptr, xray.ptr, yray.ptr, to!(int)(dim[0]),
            to!(int)(dim[1]), temp.min, temp.max);
}

void listplot3d_file(string filenamenoext, double[][] data, double[] xray = null,
        double[] yray = null)
{
    string filename=filenamenoext ~ ".png";
    if (filename.exists) remove(filename);
    auto dim = data.dims();
    if (xray is null)
        xray = arange(0., to!(double)(dim[0]), 1.);
    if (yray is null)
        yray = arange(0., to!(double)(dim[1]), 1.);
    auto temp = minmax(data);
    qplsur_file(filename, data.srg.ptr, xray.ptr, yray.ptr,
            to!(int)(dim[0]), to!(int)(dim[1]), temp.min, temp.max);
}

void listplot3d(DIdx data, double[] xray = null, double[] yray = null)
{
    auto dim = data.dims();
    if ((xray is null) && (yray is null))
    {
        xray = arange(0., to!(double)(dim[0]), 1.);
        yray = arange(0., to!(double)(dim[1]), 1.);
    }
    auto temp = minmax(data);
    qplsur(data.srg.ptr, xray.ptr, yray.ptr, to!(int)(dim[0]),
            to!(int)(dim[1]), temp.min, temp.max);
}

void listplot3d_file(string filenamenoext, DIdx data, double[] xray = null, double[] yray = null)
{
    string filename=filenamenoext ~ ".png";
    if (filename.exists) remove(filename);
    auto dim = data.dims();
    if ((xray is null) && (yray is null))
    {
        xray = arange(0., to!(double)(dim[0]), 1.);
        yray = arange(0., to!(double)(dim[1]), 1.);
    }
    auto temp = minmax(data);
    qplsur_file(filename, data.srg.ptr, xray.ptr, yray.ptr,
            to!(int)(dim[0]), to!(int)(dim[1]), temp.min, temp.max);
}

void densityplot(double[][] data, double[] xray = null, double[] yray = null)
{
    auto dim = data.dims();
    if ((xray is null) && (yray is null))
    {
        xray = arange(0., to!(double)(dim[0]), 1.);
        yray = arange(0., to!(double)(dim[1]), 1.);
    }
    auto temp = minmax(data);
    qplcrv(data.srg.ptr, xray.ptr, yray.ptr, to!(int)(dim[0]),
            to!(int)(dim[1]), temp.min, temp.max);
}

void densityplot(DIdx data, double[] xray = null, double[] yray = null)
{
    auto dim = data.dims();
    if (xray is null)
        xray = arange(0., to!(double)(dim[0]), 1.);
    if (yray is null)
        yray = arange(0., to!(double)(dim[1]), 1.);
    auto temp = minmax(data);
    qplcrv(data.srg.ptr, xray.ptr, yray.ptr, to!(int)(dim[0]),
            to!(int)(dim[1]), temp.min, temp.max);
}

void densityplot_file(string filenamenoext, double[][] data,
        double[] xray = null, double[] yray = null)
{
    string filename=filenamenoext ~ ".png";
    if (filename.exists) remove(filename);
    auto dim = data.dims();
    if ((xray is null) && (yray is null))
    {
        xray = arange(0., to!(double)(dim[0]), 1.);
        yray = arange(0., to!(double)(dim[1]), 1.);
    }
    auto temp = minmax(data);
    qplcrv_file(filename, data.srg.ptr, xray.ptr, yray.ptr,
            to!(int)(dim[0]), to!(int)(dim[1]), temp.min, temp.max);
}

void densityplot_file(string filenamenoext, DIdx data, double[] xray = null, double[] yray = null)
{
    string filename=filenamenoext ~ ".png";
    if (filename.exists) remove(filename);
    auto dim = data.dims();
    if (xray is null)
        xray = arange(0., to!(double)(dim[0]), 1.);
    if (yray is null)
        yray = arange(0., to!(double)(dim[1]), 1.);
    auto temp = minmax(data);
    qplcrv_file(filename, data.srg.ptr, xray.ptr, yray.ptr,
            to!(int)(dim[0]), to!(int)(dim[1]), temp.min, temp.max);
}

void qplsur(double* zdata, double* xray, double* yray, int nx, int ny, double zmin, double zmax)
{
    metafl(cast(char*) "XWIN");
    setpag(cast(char*) "DA4L");
    disini();
    winkey(cast(char*) "RETURN");
//    winkey(cast(char*) "ESCAPE");
    pagera();
    complx();
    axspos(400, 1800);
    axslen(2400, 1600);
    name(cast(char*) "X-axis", cast(char*) "x");
    name(cast(char*) "Y-axis", cast(char*) "y");
    name(cast(char*) "Z-axis", cast(char*) "z");
    graf3d(xray[0], xray[nx - 1], xray[0], 0.2 * (xray[nx - 1] - xray[0]),
            yray[0], yray[ny - 1], yray[0], 0.2 * (xray[nx - 1] - xray[0]), zmin,
            zmax, zmin, 0.2 * (zmax - zmin));
    height(50);
    shdmod(cast(char*) "smooth", cast(char*) "surface");
    surshd(xray, nx, yray, ny, zdata);
    disfin();
}

void qplsur_file(string filename, double* zdata, double* xray, double* yray,
        int nx, int ny, double zmin, double zmax)
{
    metafl(cast(char*) "PNG");
    setfil(cast(char*) filename);
    disini();
    winkey(cast(char*) "ESCAPE");
    pagera();
    complx();
    axspos(400, 1800);
    axslen(2400, 1600);
    name(cast(char*) "X-axis", cast(char*) "x");
    name(cast(char*) "Y-axis", cast(char*) "y");
    name(cast(char*) "Z-axis", cast(char*) "z");
    graf3d(xray[0], xray[nx - 1], xray[0], 0.2 * (xray[nx - 1] - xray[0]),
            yray[0], yray[ny - 1], yray[0], 0.2 * (xray[nx - 1] - xray[0]), zmin,
            zmax, zmin, 0.2 * (zmax - zmin));
    height(50);
    shdmod(cast(char*) "smooth", cast(char*) "surface");
    surshd(xray, nx, yray, ny, zdata);
    disfin();
}

void qplcrv(double* z, double* xray, double* yray, int nx, int ny, double zmin, double zmax)
{
    metafl(cast(char*) "XWIN");
    setpag(cast(char*) "DA4L");
    disini();
    //winkey(cast(char*) "ESCAPE");
    winkey(cast(char*) "RETURN");
    pagera();
    complx();
    axspos(400, 1800);
    ax3len(2000, 1600, 1500);
    graf3(xray[0], xray[nx - 1], xray[0], 0.2 * (xray[nx - 1] - xray[0]),
            yray[0], yray[ny - 1], yray[0], 0.2 * (xray[nx - 1] - xray[0]), zmin,
            zmax, zmin, 0.2 * (zmax - zmin));
    crvmat(z, nx, ny, 1, 1);
    disfin();
}

void qplcrv_file(string filename, double* z, double* xray, double* yray, int nx,
        int ny, double zmin, double zmax)
{
    metafl(cast(char*) "PNG");
    setfil(cast(char*) filename);
    disini();
    winkey(cast(char*) "ESCAPE");
    pagera();
    complx();
    axspos(400, 1800);
    ax3len(2000, 1600, 1500);
    graf3(xray[0], xray[nx - 1], xray[0], 0.2 * (xray[nx - 1] - xray[0]),
            yray[0], yray[ny - 1], yray[0], 0.2 * (xray[nx - 1] - xray[0]), zmin,
            zmax, zmin, 0.2 * (zmax - zmin));
    crvmat(z, nx, ny, 1, 1);
    disfin();
}

void movieplot_begin()
{
    metafl(cast(char*) "XWIN");
    setpag(cast(char*) "DA4L");
    disini();
    pagera();
//    winkey(cast(char*) "ESCAPE");
    winkey(cast(char*) "RETURN");
    complx();
}

void paperplot_begin(string filename)
{
    metafl(cast(char*) "PNG");
    setfil(cast(char*) filename);
    disini();
    pagera();
    winkey(cast(char*) "ESCAPE");
    complx();
}

void listplot_m(U)(U[] y, double[] x)
{
    static if (is(U == double))
    {
        int n = cast(int) y.length;
        double[] xg, data;
        xg = x;
        data = y;
        auto temp = minmax(y);
        double ymin, ymax;
        ymin = temp.min;
        ymax = temp.max;
        endgrf();
        erase();
        axspos(400, 1800);
        axslen(2400, 1600);
        name(cast(char*) "X-axis", cast(char*) "x");
        name(cast(char*) "Y-axis", cast(char*) "y");
        color(cast(char*) "fore");
        titlin(cast(char*) "", 1);
        graf(xg[0], xg[n - 1], xg[0], (xg[n - 1] - xg[0]) / 5., ymin, ymax,
                ymin, (ymax - ymin) / 5.);
        title();
        color(cast(char*) "white");
        curve(xg.ptr, data.ptr, cast(int) xg.length);
    }
    else static if (is(U : Complex!double))
    {
        size_t n = y.length;
        double[] xg, data, argdata;
        xg = x;
        data.length = n;
        argdata.length = n;
        for (size_t i = 0; i < n; ++i)
        {
            data[i] = cabs(y[i]);
            argdata[i] = atan2(y[i].im, y[i].re) / double(PI);
        }

        auto temp = minmax(data);
        double ymin = temp.min, ymax = temp.max;
        for (size_t i = 0; i < n; ++i)
        {
            data[i] = data[i] / ymax;
        }
        endgrf();
        erase();
        color(cast(char*) "fore");

        titlin(cast(char*) "Normalized abs (white) and arg (red)", 1);
        version (D_Version2)
            titlin(cast(char*)("normalized to absmax=" ~ to!(string)(ymax) ~ "; argmax=Pi"), 3);
        else
            titlin(cast(char*)("normalized to absmax=" ~ to!(string)(ymax) ~ "; argmax=Pi"), 3);

        graf(xg[0], xg[n - 1], xg[0], (xg[n - 1] - xg[0]) / 5., -1., 1., -1., 0.2);
        title();
        color(cast(char*) "white");
        curve(xg.ptr, data.ptr, cast(int) xg.length);
        color(cast(char*) "red");
        curve(xg.ptr, argdata.ptr, cast(int) xg.length);
    }
    else
        assert(false);
}

void listplot_m(U)(Idx!U y, double[] x)
{
  listplot_m(y.srg, x);
}

void plot_end()
{
    disfin();
}

alias plot_end movieplot_end;
