/**
 * Basic vector, matrix, idx tools.
 *
 * Copyright:   Copyright (C) 2021 by Yinchieh Lai
 * Author:     Yinchieh Lai
 * License:     $(LINK2 http://www.boost.org/LICENSE_1_0.txt, Boost License 1.0)
 */

module addon.speigen;
import addon.tool;
import addon.spm;
import addon.splinear;
import addon.umfpack : umfpack_di_free_numeric, umfpack_zi_free_numeric;


pragma(lib, "arpack");
pragma(lib, "umfpack");

extern (C)
{
    void dsaupd_(int* ido, char* bmat, int* n, char* which, int* nev,
            double* tol, double* resid, int* ncv, double* v, int* ldv,
            int* iparam, int* ipntr, double* workd, double* workl, int* lworkl, int* info);

    void dseupd_(int* rvec, char* All, int* select, double* d, double* z,
            int* ldz, double* sigma, char* bmat, int* n, char* which, int* nev,
            double* tol, double* resid, int* ncv, double* v, int* ldv,
            int* iparam, int* ipntr, double* workd, double* workl, int* lworkl, int* ierr);

    void znaupd_(int* ido, char* bmat, int* n, char* which, int* nev,
            double* tol, Cdouble* resid, int* ncv, Cdouble* v, int* ldv,
            int* iparam, int* ipntr, Cdouble* workd, Cdouble* workl,
            int* lworkl, double* rwork, int* info);

    void zneupd_(int* rvec, char* All, int* select, Cdouble* d, Cdouble* z,
            int* ldz, double* sigma, Cdouble* workev, char* bmat, int* n,
            char* which, int* nev, double* tol, Cdouble* resid, int* ncv,
            Cdouble* v, int* ldv, int* iparam, int* ipntr, Cdouble* workd,
            Cdouble* workl, int* lworkl, double* rwork, int* ierr);
}

double[] eigenvalues_ds(Spm!(double) M0, int nev, double shift, double[] Evals = null)
{
    auto M = M0.dup;
    for (size_t i = 0; i < M.dim[0]; ++i)
        M[i, i] -= shift;
    auto m = M.compact();
    int n = cast(int) m.count.length - 1;
    if (Evals is null)
        Evals = newvec(nev, double.init);
    int ido = 0; /* Initialization of the reverse communication
		  parameter. */

    char[1] bmat = "I"; /* Specifies that the right hand side matrix
			 should be the identity matrix; this makes
			 the writelnoblem a standard eigenvalue writelnoblem.
			 Setting bmat = "G" would have us solve the
			 writelnoblem Av = lBv (this would involve using
			 some other writelnograms from BLAS, however). */

    char[2] which = "LM";
    /* Ask for the nev eigenvalues of smallest
			   magnitude.  The possible options are
			   LM: largest magnitude
			   SM: smallest magnitude
			   LA: largest real component
			   SA: smallest real compoent
			   LI: largest imaginary component
			   SI: smallest imaginary component */

    double tol = 0.0; /* Sets the tolerance; tol<=0 specifies 
		       machine writelnecision */

    double[] resid;
    resid.length = n;

    int ncv = 2 * nev; /* The largest number of basis vectors that will
		      be used in the Implicitly Restarted Arnoldi
		      writelnocess.  Work per major iteration is
		      writelnoportional to N*NCV*NCV. */
    if (ncv > n)
        ncv = n;
    if (ncv < 12)
        ncv = 12;

    double[] v;
    int ldv = n;
    v.length = ldv * ncv;

    int[] iparam;
    iparam.length = 11; /* An array used to pass information to the routines
			   about their functional modes. */
    iparam[0] = 1; // Specifies the shift strategy (1->exact)
    iparam[2] = 3 * n; // Maximum number of iterations
    iparam[6] = 1; /* Sets the mode of dsaupd.
		      1 is exact shifting,
		      2 is user-supplied shifts,
		      3 is shift-invert mode,
		      4 is buckling mode,
		      5 is Cayley mode. */

    int[] ipntr;
    ipntr.length = 11; /* Indicates the locations in the work array workd
			  where the input and output vectors in the
			  callback routine are located. */

    double[] workd;
    workd.length = 3 * n;

    double[] workl;
    workl.length = ncv * (ncv + 8);

    int lworkl = ncv * (ncv + 8); /* Length of the workl array */

    int info = 0; /* Passes convergence information out of the iteration
		   routine. */

    int rvec = 0; /* Specifies that eigenvectors should not be calculated */

    int[] select;
    select.length = ncv;
    double[] d;
    d.length = 2 * ncv; /* This vector will return the eigenvalues from
			    the second routine, dseupd. */
    double sigma;
    int ierr;
    double[] vin, vout;

    int uflag = 0;
    void* numeric;

    /* Here we enter the main loop where the calculations are
     performed.  The communication parameter ido tells us when
     the desired tolerance is reached, and at that point we exit
     and extract the solutions. */

    do
    {
        dsaupd_(&ido, bmat.ptr, &n, which.ptr, &nev, &tol, resid.ptr, &ncv,
                v.ptr, &ldv, iparam.ptr, ipntr.ptr, workd.ptr, workl.ptr, &lworkl, &info);
        if ((ido == 1) || (ido == -1))
        {
            vin = (workd.ptr + ipntr[0] - 1)[0 .. n];
            vout = (workd.ptr + ipntr[1] - 1)[0 .. n];
            if (uflag == 0)
            {
                linearsolve_store_numeric(m, numeric, vin, vout);
                uflag = 1;
            }
            else
                linearsolve_with_numeric(m, numeric, vin, vout);
        }
    }
    while ((ido == 1) || (ido == -1));

    /* From those results, the eigenvalues and vectors are
     extracted. */

    if (info < 0)
    {
        writeln("Error with dsaupd, info = ", info);
        writeln("Check documentation in dsaupd");
    }
    else
    {
        dseupd_(&rvec, cast(char*) "All".ptr, select.ptr, d.ptr, v.ptr,
                &ldv, &sigma, bmat.ptr, &n, which.ptr, &nev, &tol,
                resid.ptr, &ncv, v.ptr, &ldv, iparam.ptr, ipntr.ptr,
                workd.ptr, workl.ptr, &lworkl, &ierr);

        if (ierr != 0)
        {
            writeln("Error with dseupd, info = ", ierr);
            writeln("Check the documentation of dseupd.");
        }
        else if (info == 1)
        {
            writeln("Maximum number of iterations reached.");
        }
        else if (info == 3)
        {
            writeln("No shifts could be applied during implicit");
            writeln("Arnoldi update, try increasing NCV.");
        }

        /* Before exiting, we copy the solution information over to
       the arrays of the calling writelnogram, then clean up the
       memory used by this routine.  For some reason, when I
       don't find the eigenvectors I need to reverse the order of
       the values. */

        for (size_t i = 0; i < nev; i++)
            Evals[i] = 1. / d[nev - 1 - i] + shift;
    }
    umfpack_di_free_numeric(&numeric);
    resid.destroy();
    v.destroy();
    iparam.destroy();
    ipntr.destroy();
    workd.destroy();
    workl.destroy();
    select.destroy();
    d.destroy();

    return Evals;
}

Cdouble[] eigenvalues(Spm!(Cdouble) M0, int nev, Cdouble shift, Cdouble[] Evals = null)
{
    auto M = M0.dup;
    for (size_t i = 0; i < M.dim[0]; ++i)
        M[i, i] -= shift;
    auto m = M.compact_ri();
    int n = cast(int) m.count.length - 1;
    if (Evals is null)
        Evals = newvec(nev, Cdouble.init);
    int ido = 0; /* Initialization of the reverse communication
		  parameter. */

    char[1] bmat = "I"; /* Specifies that the right hand side matrix
			 should be the identity matrix; this makes
			 the writelnoblem a standard eigenvalue writelnoblem.
			 Setting bmat = "G" would have us solve the
			 writelnoblem Av = lBv (this would involve using
			 some other writelnograms from BLAS, however). */

    char[2] which = "LM"; /* Ask for the nev eigenvalues of smallest
			   magnitude.  The possible options are
			   LM: largest magnitude
			   SM: smallest magnitude
			   LA: largest real component
			   SA: smallest real compoent
			   LI: largest imaginary component
			   SI: smallest imaginary component */

    double tol = 0.0; /* Sets the tolerance; tol<=0 specifies 
		       machine writelnecision */

    Cdouble[] resid;
    resid.length = n;

    int ncv = 4 * nev; /* The largest number of basis vectors that will
		      be used in the Implicitly Restarted Arnoldi
		      writelnocess.  Work per major iteration is
		      writelnoportional to N*NCV*NCV. */
    if (ncv > n)
        ncv = n;
    if (ncv < 12)
        ncv = 12;

    Cdouble[] v;
    int ldv = n;
    v.length = ldv * ncv;

    int[] iparam;
    iparam.length = 11; /* An array used to pass information to the routines
			   about their functional modes. */
    iparam[0] = 1; // Specifies the shift strategy (1->exact)
    iparam[2] = 3 * n; // Maximum number of iterations
    iparam[6] = 1; /* Sets the mode of dsaupd.
		      1 is exact shifting,
		      2 is user-supplied shifts,
		      3 is shift-invert mode,
		      4 is buckling mode,
		      5 is Cayley mode. */

    int[] ipntr;
    ipntr.length = 14; /* Indicates the locations in the work array workd
			  where the input and output vectors in the
			  callback routine are located. */

    Cdouble[] workd;
    workd.length = 3 * n;

    int lworkl = 3 * ncv * ncv + 5 * ncv; /* Length of the workl array */
    Cdouble[] workl;
    workl.length = lworkl;

    double[] rwork;
    rwork.length = ncv;

    int info = 0; /* Passes convergence information out of the iteration
		   routine. */

    int rvec = 0; /* Specifies that eigenvectors should not be calculated */

    int[] select;
    select.length = ncv;
    Cdouble[] d;
    d.length = ncv; /* This vector will return the
				     eigenvalues from the second routine,
				     dseupd. */
    double sigma;

    Cdouble[] workev;
    workev.length = 3 * ncv; // I don't know what this is used for

    int ierr;

    int uflag = 0;
    void* numeric;
    /* Here we enter the main loop where the calculations are
     performed.  The communication parameter ido tells us when
     the desired tolerance is reached, and at that point we exit
     and extract the solutions. */
    double[] cvin, cvout;
    double[] vin, vout;
    vin.length = 2 * n;
    vout.length = 2 * n;
    size_t j, k;

    do
    {
        znaupd_(&ido, bmat.ptr, &n, which.ptr, &nev, &tol, resid.ptr, &ncv,
                v.ptr, &ldv, iparam.ptr, ipntr.ptr, workd.ptr, workl.ptr,
                &lworkl, rwork.ptr, &info);

        if ((ido == 1) || (ido == -1))
        {
            cvin = (cast(double*)(workd.ptr + ipntr[0] - 1))[0 .. 2 * n];
            cvout = (cast(double*)(workd.ptr + ipntr[1] - 1))[0 .. 2 * n];
            j = 0;
            k = 1;
            for (size_t i = 0; i < n; ++i)
            {
                vin[i] = cvin[j];
                vin[i + n] = cvin[k];
                j += 2;
                k += 2;
            }
            if (uflag == 0)
            {
                linearsolve_store_numeric(m, numeric, vin, vout);
                uflag = 1;
            }
            else
                linearsolve_with_numeric(m, numeric, vin, vout);
            j = 0;
            k = 1;
            for (size_t i = 0; i < n; ++i)
            {
                cvout[j] = vout[i];
                cvout[k] = vout[i + n];
                j += 2;
                k += 2;
            }
        }
    }
    while ((ido == 1) || (ido == -1));

    /* From those results, the eigenvalues and vectors are
     extracted. */

    if (info < 0)
    {
        writeln("Error with znaupd, info = ", info);
        writeln("Check documentation in dsaupd");
    }
    else
    {
        zneupd_(&rvec, cast(char*) "All".ptr, select.ptr, d.ptr, v.ptr, &ldv,
                &sigma, workev.ptr, bmat.ptr, &n, which.ptr, &nev, &tol,
                resid.ptr, &ncv, v.ptr, &ldv, iparam.ptr, ipntr.ptr,
                workd.ptr, workl.ptr, &lworkl, rwork.ptr, &ierr);

        if (ierr != 0)
        {
            writeln("Error with zneupd, info = ", ierr);
            writeln("Check the documentation of zneupd.");
        }
        else if (info == 1)
        {
            writeln("Maximum number of iterations reached.");
        }
        else if (info == 3)
        {
            writeln("No shifts could be applied during implicit");
            writeln("Arnoldi update, try increasing NCV.");
        }

        /* Before exiting, we copy the solution information over to
       the arrays of the calling writelnogram, then clean up the
       memory used by this routine.  For some reason, when I
       don't find the eigenvectors I need to reverse the order of
       the values. */

        for (size_t i = 0; i < nev; i++)
            Evals[i] = complex(1.0, 0.0) / d[nev - 1 - i] + shift;
    }
    // Sort the energies by real part
    umfpack_zi_free_numeric(&numeric);
    resid.destroy();
    v.destroy();
    iparam.destroy();
    ipntr.destroy();
    workd.destroy();
    workl.destroy();
    select.destroy();
    d.destroy();
    vin.destroy();
    vout.destroy();

    return Evals;
}

Tuple!(double[], "eigenvalues", double[][], "eigenvectors") eigensystem_ds(
        Spm!(double) M0, int nev, double shift, double[] Evals = null, double[][] Evecs = null)
{
    auto M = M0.dup;
    for (size_t i = 0; i < M.dim[0]; ++i)
        M[i, i] -= shift;
    auto m = M.compact();
    int n = cast(int) m.count.length - 1;
    if (Evals is null)
        Evals = newvec(nev, double.init);
    if (Evecs is null)
        Evecs = newarray!(double)(nev, n);
    int ido = 0; /* Initialization of the reverse communication
		  parameter. */

    char[1] bmat = "I"; /* Specifies that the right hand side matrix
			 should be the identity matrix; this makes
			 the writelnoblem a standard eigenvalue writelnoblem.
			 Setting bmat = "G" would have us solve the
			 writelnoblem Av = lBv (this would involve using
			 some other writelnograms from BLAS, however). */

    char[2] which = "LM";
    /* Ask for the nev eigenvalues of smallest
			   magnitude.  The possible options are
			   LM: largest magnitude
			   SM: smallest magnitude
			   LA: largest real component
			   SA: smallest real compoent
			   LI: largest imaginary component
			   SI: smallest imaginary component */

    double tol = 0.0; /* Sets the tolerance; tol<=0 specifies 
		       machine writelnecision */

    double[] resid;
    resid.length = n;

    int ncv = 4 * nev; /* The largest number of basis vectors that will
		      be used in the Implicitly Restarted Arnoldi
		      writelnocess.  Work per major iteration is
		      writelnoportional to N*NCV*NCV. */
    if (ncv > n)
        ncv = n;
    if (ncv < 12)
        ncv = 12;

    double[] v;
    int ldv = n;
    v.length = ldv * ncv;

    int[] iparam;
    iparam.length = 11; /* An array used to pass information to the routines
			   about their functional modes. */
    iparam[0] = 1; // Specifies the shift strategy (1->exact)
    iparam[2] = 3 * n; // Maximum number of iterations
    iparam[6] = 1; /* Sets the mode of dsaupd.
		      1 is exact shifting,
		      2 is user-supplied shifts,
		      3 is shift-invert mode,
		      4 is buckling mode,
		      5 is Cayley mode. */

    int[] ipntr;
    ipntr.length = 11; /* Indicates the locations in the work array workd
			  where the input and output vectors in the
			  callback routine are located. */

    double[] workd;
    workd.length = 3 * n;

    double[] workl;
    workl.length = ncv * (ncv + 8);

    int lworkl = ncv * (ncv + 8); /* Length of the workl array */

    int info = 0; /* Passes convergence information out of the iteration
		   routine. */

    int rvec = 1; /* Specifies that eigenvectors should not be calculated */

    int[] select;
    select.length = ncv;
    double[] d;
    d.length = 2 * ncv; /* This vector will return the eigenvalues from
			    the second routine, dseupd. */
    double sigma;
    int ierr;

    int uflag = 0;
    void* numeric;
    /* Here we enter the main loop where the calculations are
     performed.  The communication parameter ido tells us when
     the desired tolerance is reached, and at that point we exit
     and extract the solutions. */

    double[] vin, vout;
    do
    {
        dsaupd_(&ido, bmat.ptr, &n, which.ptr, &nev, &tol, resid.ptr, &ncv,
                v.ptr, &ldv, iparam.ptr, ipntr.ptr, workd.ptr, workl.ptr, &lworkl, &info);

        if ((ido == 1) || (ido == -1))
        {
            vin = (workd.ptr + ipntr[0] - 1)[0 .. n];
            vout = (workd.ptr + ipntr[1] - 1)[0 .. n];
            if (uflag == 0)
            {
                linearsolve_store_numeric(m, numeric, vin, vout);
                uflag = 1;
            }
            else
                linearsolve_with_numeric(m, numeric, vin, vout);
        }
    }
    while ((ido == 1) || (ido == -1));

    /* From those results, the eigenvalues and vectors are
     extracted. */

    if (info < 0)
    {
        writeln("Error with dsaupd, info = ", info);
        writeln("Check documentation in dsaupd");
    }
    else
    {
        dseupd_(&rvec, cast(char*) "All".ptr, select.ptr, d.ptr, v.ptr,
                &ldv, &sigma, bmat.ptr, &n, which.ptr, &nev, &tol,
                resid.ptr, &ncv, v.ptr, &ldv, iparam.ptr, ipntr.ptr,
                workd.ptr, workl.ptr, &lworkl, &ierr);

        if (ierr != 0)
        {
            writeln("Error with dseupd, info = ", ierr);
            writeln("Check the documentation of dseupd.");
        }
        else if (info == 1)
        {
            writeln("Maximum number of iterations reached.");
        }
        else if (info == 3)
        {
            writeln("No shifts could be applied during implicit");
            writeln("Arnoldi update, try increasing NCV.");
        }

        /* Before exiting, we copy the solution information over to
       the arrays of the calling writelnogram, then clean up the
       memory used by this routine.  For some reason, when I
       don't find the eigenvectors I need to reverse the order of
       the values. */

        /*    for (int i=0; i<nev; i++) Evals[i] = d[i];
    for (int i=0; i<nev; i++) for (int j=0; j<n; j++) Evecs[i,j] = v[i*n+j]; */

        for (int i = 0; i < nev; i++)
            Evals[i] = 1. / d[i] + shift;
        Evecs.srg[] = v[0 .. n * nev];

    }
    Tuple!(double[], "eigenvalues", double[][], "eigenvectors") result;
    result.eigenvalues = Evals;
    result.eigenvectors = Evecs;

    umfpack_di_free_numeric(&numeric);
    resid.destroy();
    v.destroy();
    iparam.destroy();
    ipntr.destroy();
    workd.destroy();
    workl.destroy();
    select.destroy();
    d.destroy();

    return result;
}

Tuple!(Cdouble[], "eigenvalues", Cdouble[][], "eigenvectors") eigensystem(Spm!(
        Cdouble) M0, int nev, Cdouble shift, Cdouble[] Evals = null, Cdouble[][] Evecs = null)
{
    auto M = M0.dup;
    for (size_t i = 0; i < M.dim[0]; ++i)
        M[i, i] -= shift;
    auto m = M.compact_ri();
    int n = cast(int) m.count.length - 1;
    if (Evals is null)
        Evals = newvec(nev, Cdouble.init);
    if (Evecs is null)
        Evecs = newarray!(Cdouble)(nev, n);
    int ido = 0;
    char[1] bmat = "I";
    char[2] which = "LM";
    double tol = 0.0;
    Cdouble[] resid;
    resid.length = n;
    int ncv = 4 * nev;
    if (ncv > n)
        ncv = n;
    if (ncv < 12)
        ncv = 12;
    Cdouble[] v;
    int ldv = n;
    v.length = ldv * ncv;
    int[] iparam;
    iparam.length = 11;
    iparam[0] = 1;
    iparam[2] = 3 * n;
    iparam[6] = 1;
    int[] ipntr;
    ipntr.length = 14;
    Cdouble[] workd;
    workd.length = 3 * n;
    int lworkl = 3 * ncv * ncv + 5 * ncv;
    Cdouble[] workl;
    workl.length = lworkl;
    double[] rwork;
    rwork.length = ncv;
    int info = 0;
    int rvec = 1; // Changed from above
    int[] select;
    select.length = ncv;
    Cdouble[] d;
    d.length = ncv;
    double sigma;
    Cdouble[] workev;
    workev.length = 3 * ncv;
    int ierr;

    int uflag = 0;
    void* numeric;

    double[] cvin, cvout;
    double[] vin, vout;
    vin.length = 2 * n;
    vout.length = 2 * n;
    size_t j, k;

    do
    {
        znaupd_(&ido, bmat.ptr, &n, which.ptr, &nev, &tol, resid.ptr, &ncv,
                v.ptr, &ldv, iparam.ptr, ipntr.ptr, workd.ptr, workl.ptr,
                &lworkl, rwork.ptr, &info);

        if ((ido == 1) || (ido == -1))
        {
            cvin = (cast(double*)(workd.ptr + ipntr[0] - 1))[0 .. 2 * n];
            cvout = (cast(double*)(workd.ptr + ipntr[1] - 1))[0 .. 2 * n];
            j = 0;
            k = 1;
            for (size_t i = 0; i < n; ++i)
            {
                vin[i] = cvin[j];
                vin[i + n] = cvin[k];
                j += 2;
                k += 2;
            }
            if (uflag == 0)
            {
                linearsolve_store_numeric(m, numeric, vin, vout);
                uflag = 1;
            }
            else
                linearsolve_with_numeric(m, numeric, vin, vout);
            j = 0;
            k = 1;
            for (size_t i = 0; i < n; ++i)
            {
                cvout[j] = vout[i];
                cvout[k] = vout[i + n];
                j += 2;
                k += 2;
            }
        }
    }
    while ((ido == 1) || (ido == -1));

    if (info < 0)
    {
        writeln("Error with znaupd, info = ", info);
        writeln("Check documentation in dsaupd");
    }
    else
    {
        zneupd_(&rvec, cast(char*) "All".ptr, select.ptr, d.ptr, v.ptr, &ldv,
                &sigma, workev.ptr, bmat.ptr, &n, which.ptr, &nev, &tol,
                resid.ptr, &ncv, v.ptr, &ldv, iparam.ptr, ipntr.ptr,
                workd.ptr, workl.ptr, &lworkl, rwork.ptr, &ierr);

        if (ierr != 0)
        {
            writeln("Error with zneupd, info = ", ierr);
            writeln("Check the documentation of zneupd.");
        }
        else if (info == 1)
        {
            writeln("Maximum number of iterations reached.");
        }
        else if (info == 3)
        {
            writeln("No shifts could be applied during implicit");
            writeln("Arnoldi update, try increasing NCV.");
        }
        //    
        //    size_t i, j;
        //    for (i=0; i<nev; i++) Evals[i] = d[i];
        //    for (i=0; i<nev; i++) for (j=0; j<n; j++) Evecs[i,j] = v[i*n+j]; 
        //    

        //    Evals[]=d[0..nev];
        for (size_t i = 0; i < nev; i++)
            Evals[i] = 1. / d[i] + shift;
        Evecs.srg[] = v[0 .. n * nev];
        /*
    Cdouble temp;
    for (i=0; i<nev; i++) {
      for (j=i; j<nev; j++) {
	if (Evals[j].re > Evals[i].re) {
	  temp = Evals[j];
	  Evals[j] = Evals[i];
	  Evals[i] = temp;
	  for (int k=0; k<n; k++) {
	    temp = Evecs[k][i];
	    Evecs[k][i] = Evecs[k][j];
	    Evecs[k][j] = temp;
	  }
	}
      }
    }
*/
    }
    Tuple!(Cdouble[], "eigenvalues", Cdouble[][], "eigenvectors") result;
    result.eigenvalues = Evals;
    result.eigenvectors = Evecs;

    umfpack_zi_free_numeric(&numeric);
    resid.destroy();
    v.destroy();
    iparam.destroy();
    ipntr.destroy();
    workd.destroy();
    workl.destroy();
    select.destroy();
    d.destroy();
    vin.destroy();
    vout.destroy();

    return result;
}
