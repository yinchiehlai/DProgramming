/*==========================================================================
 * umfpack.d
 *    Written in the D Programming Language (http://www.digitalmars.com/d)
 */
/***************************************************************************
 * Interface to the UMFPACK sparse matrix library.
 *
 * This is a port of the UMFPACK 5.1 header file to D.
 *
 * To obtain UMFPACK itself see:
 *    http://www.cise.ufl.edu/research/sparse/umfpack
 *
 * Authors:  William V. Baxter III
 * Date: 23 Feb 2008
 * Copyright: (C) 2007-2008 William Baxter, OLM Digital, Inc.
 * License: ZLIB/LIBPNG
 */
//===========================================================================

/* -------------------------------------------------------------------------- */
/* UMFPACK is Copyright (c) Timothy A. Davis, CISE,                           */
/* Univ. of Florida.  All Rights Reserved.  See ../Doc/License for License.   */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

module addon.umfpack;

// Set this version if you compiled UMFPACK with 64-bit longs
// 'long' seems to be the default, which is 32-bits on 32-bit machines.
version (USING_UF_LONG64)
{
    alias long UF_long;
}
else
{
    alias int UF_long;
}

version (build)
{
    pragma(link, "umfpack");
    version (Linux) pragma(link, "amd");
    pragma(link, "blas");
}

enum
{
    UMFPACK_INFO = 90,
    UMFPACK_CONTROL = 20
}

extern (C)
{
    /* ========================================================================== */
    /* === umfpack_symbolic ===================================================== */
    /* ========================================================================== */

    int umfpack_di_symbolic(int n_row, int n_col,/*const*/
            int* Ap,/*const*/
            int* Ai,/*const*/
            double* Ax,
            void** Symbolic,/*const*/
            double* Control, double* Info);

    UF_long umfpack_dl_symbolic(UF_long n_row, UF_long n_col,/*const*/
            UF_long* Ap,
            /*const*/
            UF_long* Ai,/*const*/
            double* Ax, void** Symbolic,/*const*/
            double* Control, double* Info);

    int umfpack_zi_symbolic(int n_row, int n_col,/*const*/
            int* Ap,/*const*/
            int* Ai,/*const*/
            double* Ax, /*const*/ double* Az, void** Symbolic,/*const*/
            double* Control, double* Info);

    UF_long umfpack_zl_symbolic(UF_long n_row, UF_long n_col,/*const*/
            UF_long* Ap,
            /*const*/
            UF_long* Ai,/*const*/
            double* Ax, /*const*/ double* Az, void** Symbolic,/*const*/
            double* Control, double* Info,);

    /* ========================================================================== */
    /* === umfpack_numeric ====================================================== */
    /* ========================================================================== */

    int umfpack_di_numeric(/*const*/
            int* Ap,/*const*/
            int* Ai,/*const*/
            double* Ax, void* Symbolic,
            void** Numeric,/*const*/
            double* Control, double* Info,);

    UF_long umfpack_dl_numeric(/*const*/
            UF_long* Ap,/*const*/
            UF_long* Ai,/*const*/
            double* Ax,
            void* Symbolic, void** Numeric,/*const*/
            double* Control, double* Info,);

    int umfpack_zi_numeric(/*const*/
            int* Ap,/*const*/
            int* Ai,/*const*/
            double* Ax, /*const*/ double* Az,
            void* Symbolic, void** Numeric,/*const*/
            double* Control, double* Info);

    UF_long umfpack_zl_numeric(/*const*/
            UF_long* Ap,/*const*/
            UF_long* Ai,/*const*/
            double* Ax, /*const*/ double* Az,
            void* Symbolic, void** Numeric,/*const*/
            double* Control, double* Info);

    /* ========================================================================== */
    /* === umfpack_solve ======================================================== */
    /* ========================================================================== */

    int umfpack_di_solve(int sys,/*const*/
            int* Ap,/*const*/
            int* Ai,/*const*/
            double* Ax, double* X,
            /*const*/
            double* B, void* Numeric,/*const*/
            double* Control, double* Info);

    UF_long umfpack_dl_solve(UF_long sys,/*const*/
            UF_long* Ap,/*const*/
            UF_long* Ai,/*const*/
            double* Ax,
            double* X,/*const*/
            double* B, void* Numeric,/*const*/
            double* Control, double* Info);

    int umfpack_zi_solve(int sys,/*const*/
            int* Ap,/*const*/
            int* Ai,/*const*/
            double* Ax, /*const*/ double* Az,
            double* Xx, double* Xz,/*const*/
            double* Bx, /*const*/ double* Bz, void* Numeric,
            /*const*/
            double* Control, double* Info);

    UF_long umfpack_zl_solve(UF_long sys,/*const*/
            UF_long* Ap,/*const*/
            UF_long* Ai,/*const*/
            double* Ax, /*const*/ double* Az, double* Xx, double* Xz,/*const*/
            double* Bx, /*const*/ double* Bz, void* Numeric,/*const*/
            double* Control, double* Info);

    /* ========================================================================== */
    /* === umfpack_free_symbolic ================================================ */
    /* ========================================================================== */

    void umfpack_di_free_symbolic(void** Symbolic);

    void umfpack_dl_free_symbolic(void** Symbolic);

    void umfpack_zi_free_symbolic(void** Symbolic);

    void umfpack_zl_free_symbolic(void** Symbolic);

    /* ========================================================================== */
    /* === umfpack_free_numeric ================================================= */
    /* ========================================================================== */

    void umfpack_di_free_numeric(void** Numeric);

    void umfpack_dl_free_numeric(void** Numeric);

    void umfpack_zi_free_numeric(void** Numeric);

    void umfpack_zl_free_numeric(void** Numeric);

    /* ========================================================================== */
    /* === umfpack_defaults ===================================================== */
    /* ========================================================================== */

    void umfpack_di_defaults(double* Control /*[UMFPACK_CONTROL]*/
    );

    void umfpack_dl_defaults(double* Control /*[UMFPACK_CONTROL]*/
    );

    void umfpack_zi_defaults(double* Control /*[UMFPACK_CONTROL]*/
    );

    void umfpack_zl_defaults(double* Control /*[UMFPACK_CONTROL]*/
    );

    /* ========================================================================== */
    /* === umfpack_qsymbolic ==================================================== */
    /* ========================================================================== */

    int umfpack_di_qsymbolic(int n_row, int n_col,/*const*/
            int* Ap,/*const*/
            int* Ai,/*const*/
            double* Ax,
            /*const*/
            int* Qinit, void** Symbolic,/*const*/
            double* Control /*[UMFPACK_CONTROL]*/ , double* Info /*[UMFPACK_INFO]*/
    );

    UF_long umfpack_dl_qsymbolic(UF_long n_row, UF_long n_col,/*const*/
            UF_long* Ap,
            /*const*/
            UF_long* Ai,/*const*/
            double* Ax,/*const*/
            UF_long* Qinit, void** Symbolic,/*const*/
            double* Control /*[UMFPACK_CONTROL]*/ , double* Info /*[UMFPACK_INFO]*/
    );

    int umfpack_zi_qsymbolic(int n_row, int n_col,/*const*/
            int* Ap,/*const*/
            int* Ai,/*const*/
            double* Ax, /*const*/ double* Az,/*const*/
            int* Qinit, void** Symbolic,
            /*const*/
            double* Control /*[UMFPACK_CONTROL]*/ , double* Info /*[UMFPACK_INFO]*/
    );

    UF_long umfpack_zl_qsymbolic(UF_long n_row, UF_long n_col,/*const*/
            UF_long* Ap,
            /*const*/
            UF_long* Ai,/*const*/
            double* Ax, /*const*/ double* Az,/*const*/
            UF_long* Qinit, void** Symbolic,
            /*const*/
            double* Control /*[UMFPACK_CONTROL]*/ , double* Info /*[UMFPACK_INFO]*/
    );

    /* ========================================================================== */
    /* === umfpack_wsolve ======================================================= */
    /* ========================================================================== */

    int umfpack_di_wsolve(int sys,/*const*/
            int* Ap,/*const*/
            int* Ai,/*const*/
            double* Ax, double* X,
            /*const*/
            double* B, void* Numeric,/*const*/
            double* Control /*[UMFPACK_CONTROL]*/ , double* Info /*[UMFPACK_INFO]*/ ,
            int* Wi, double* W);

    UF_long umfpack_dl_wsolve(UF_long sys,/*const*/
            UF_long* Ap,/*const*/
            UF_long* Ai,/*const*/
            double* Ax,
            double* X,/*const*/
            double* B, void* Numeric,/*const*/
            double* Control /*[UMFPACK_CONTROL]*/ , double* Info /*[UMFPACK_INFO]*/ ,
            UF_long* Wi, double* W);

    int umfpack_zi_wsolve(int sys,/*const*/
            int* Ap,/*const*/
            int* Ai,/*const*/
            double* Ax, /*const*/ double* Az,
            double* Xx, double* Xz,/*const*/
            double* Bx, /*const*/ double* Bz, void* Numeric,
            /*const*/
            double* Control /*[UMFPACK_CONTROL]*/ , double* Info /*[UMFPACK_INFO]*/ ,
            int* Wi, double* W);

    UF_long umfpack_zl_wsolve(UF_long sys,/*const*/
            UF_long* Ap,/*const*/
            UF_long* Ai,/*const*/
            double* Ax, /*const*/ double* Az, double* Xx, double* Xz,
            /*const*/
            double* Bx, /*const*/ double* Bz, void* Numeric,/*const*/
            double* Control /*[UMFPACK_CONTROL]*/ , double* Info /*[UMFPACK_INFO]*/ ,
            UF_long* Wi, double* W);

    /* ========================================================================== */
    /* === umfpack_triplet_to_col =============================================== */
    /* ========================================================================== */

    int umfpack_di_triplet_to_col(int n_row, int n_col, int nz,/*const*/
            int* Ti,/*const*/
            int* Tj,
            /*const*/
            double* Tx, int* Ap, int* Ai, double* Ax, int* Map);

    UF_long umfpack_dl_triplet_to_col(UF_long n_row, UF_long n_col, UF_long nz,
            /*const*/
            UF_long* Ti,/*const*/
            UF_long* Tj,/*const*/
            double* Tx, UF_long* Ap, UF_long* Ai, double* Ax, UF_long* Map);

    int umfpack_zi_triplet_to_col(int n_row, int n_col, int nz,/*const*/
            int* Ti,/*const*/
            int* Tj,
            /*const*/
            double* Tx, /*const*/ double* Tz, int* Ap, int* Ai, double* Ax, double* Az, int* Map);

    UF_long umfpack_zl_triplet_to_col(UF_long n_row, UF_long n_col, UF_long nz,
            /*const*/
            UF_long* Ti,/*const*/
            UF_long* Tj,/*const*/
            double* Tx, /*const*/ double* Tz, UF_long* Ap,
            UF_long* Ai, double* Ax, double* Az, UF_long* Map);

    /* ========================================================================== */
    /* === umfpack_col_to_triplet =============================================== */
    /* ========================================================================== */

    int umfpack_di_col_to_triplet(int n_col,/*const*/
            int* Ap, int* Tj);

    UF_long umfpack_dl_col_to_triplet(UF_long n_col,/*const*/
            UF_long* Ap, UF_long* Tj);

    int umfpack_zi_col_to_triplet(int n_col,/*const*/
            int* Ap, int* Tj);

    UF_long umfpack_zl_col_to_triplet(UF_long n_col,/*const*/
            UF_long* Ap, UF_long* Tj);

    /* ========================================================================== */
    /* === umfpack_transpose ==================================================== */
    /* ========================================================================== */

    int umfpack_di_transpose(int n_row, int n_col,/*const*/
            int* Ap,/*const*/
            int* Ai,/*const*/
            double* Ax,
            /*const*/
            int* P,/*const*/
            int* Q, int* Rp, int* Ri, double* Rx);

    UF_long umfpack_dl_transpose(UF_long n_row, UF_long n_col,/*const*/
            UF_long* Ap,
            /*const*/
            UF_long* Ai,/*const*/
            double* Ax,/*const*/
            UF_long* P,/*const*/
            UF_long* Q, UF_long* Rp, UF_long* Ri, double* Rx);

    int umfpack_zi_transpose(int n_row, int n_col,/*const*/
            int* Ap,/*const*/
            int* Ai,/*const*/
            double* Ax, /*const*/ double* Az,/*const*/
            int* P,/*const*/
            int* Q, int* Rp, int* Ri,
            double* Rx, double* Rz, int do_conjugate);

    UF_long umfpack_zl_transpose(UF_long n_row, UF_long n_col,/*const*/
            UF_long* Ap,
            /*const*/
            UF_long* Ai,/*const*/
            double* Ax, /*const*/ double* Az,/*const*/
            UF_long* P,/*const*/
            UF_long* Q,
            UF_long* Rp, UF_long* Ri, double* Rx, double* Rz, UF_long do_conjugate);

    /* ========================================================================== */
    /* === umfpack_scale ======================================================== */
    /* ========================================================================== */

    int umfpack_di_scale(double* X,/*const*/
            double* B, void* Numeric);

    UF_long umfpack_dl_scale(double* X,/*const*/
            double* B, void* Numeric);

    int umfpack_zi_scale(double* Xx, double* Xz,/*const*/
            double* Bx, /*const*/ double* Bz, void* Numeric);

    UF_long umfpack_zl_scale(double* Xx, double* Xz,/*const*/
            double* Bx, /*const*/ double* Bz, void* Numeric);

    /* ========================================================================== */
    /* === umfpack_get_lunz ===================================================== */
    /* ========================================================================== */

    int umfpack_di_get_lunz(int* lnz, int* unz, int* n_row, int* n_col,
            int* nz_udiag, void* Numeric);

    UF_long umfpack_dl_get_lunz(UF_long* lnz, UF_long* unz, UF_long* n_row,
            UF_long* n_col, UF_long* nz_udiag, void* Numeric);

    int umfpack_zi_get_lunz(int* lnz, int* unz, int* n_row, int* n_col,
            int* nz_udiag, void* Numeric);

    UF_long umfpack_zl_get_lunz(UF_long* lnz, UF_long* unz, UF_long* n_row,
            UF_long* n_col, UF_long* nz_udiag, void* Numeric);

    /* ========================================================================== */
    /* === umfpack_get_numeric ================================================== */
    /* ========================================================================== */

    int umfpack_di_get_numeric(int* Lp, int* Lj, double* Lx, int* Up, int* Ui,
            double* Ux, int* P, int* Q, double* Dx, int* do_recip, double* Rs, void* Numeric);

    UF_long umfpack_dl_get_numeric(UF_long* Lp, UF_long* Lj, double* Lx,
            UF_long* Up, UF_long* Ui, double* Ux, UF_long* P, UF_long* Q,
            double* Dx, UF_long* do_recip, double* Rs, void* Numeric);

    int umfpack_zi_get_numeric(int* Lp, int* Lj, double* Lx, double* Lz,
            int* Up, int* Ui, double* Ux, double* Uz, int* P, int* Q,
            double* Dx, double* Dz, int* do_recip, double* Rs, void* Numeric);

    UF_long umfpack_zl_get_numeric(UF_long* Lp, UF_long* Lj, double* Lx,
            double* Lz, UF_long* Up, UF_long* Ui, double* Ux, double* Uz,
            UF_long* P, UF_long* Q, double* Dx, double* Dz, UF_long* do_recip,
            double* Rs, void* Numeric);
    /* ========================================================================== */
    /* === umfpack_get_symbolic ================================================= */
    /* ========================================================================== */

    int umfpack_di_get_symbolic(int* n_row, int* n_col, int* n1, int* nz,
            int* nfr, int* nchains, int* P, int* Q, int* Front_npivcol, int* Front_parent,
            int* Front_1strow, int* Front_leftmostdesc, int* Chain_start,
            int* Chain_maxrows, int* Chain_maxcols, void* Symbolic);

    UF_long umfpack_dl_get_symbolic(UF_long* n_row, UF_long* n_col, UF_long* n1,
            UF_long* nz, UF_long* nfr, UF_long* nchains, UF_long* P, UF_long* Q, UF_long* Front_npivcol,
            UF_long* Front_parent, UF_long* Front_1strow, UF_long* Front_leftmostdesc,
            UF_long* Chain_start, UF_long* Chain_maxrows, UF_long* Chain_maxcols, void* Symbolic);

    int umfpack_zi_get_symbolic(int* n_row, int* n_col, int* n1, int* nz,
            int* nfr, int* nchains, int* P, int* Q, int* Front_npivcol, int* Front_parent,
            int* Front_1strow, int* Front_leftmostdesc, int* Chain_start,
            int* Chain_maxrows, int* Chain_maxcols, void* Symbolic);

    UF_long umfpack_zl_get_symbolic(UF_long* n_row, UF_long* n_col, UF_long* n1,
            UF_long* nz, UF_long* nfr, UF_long* nchains, UF_long* P, UF_long* Q, UF_long* Front_npivcol,
            UF_long* Front_parent, UF_long* Front_1strow, UF_long* Front_leftmostdesc,
            UF_long* Chain_start, UF_long* Chain_maxrows, UF_long* Chain_maxcols, void* Symbolic);

    /* ========================================================================== */
    /* === umfpack_save_numeric ================================================= */
    /* ========================================================================== */

    int umfpack_di_save_numeric(void* Numeric, char* filename);

    UF_long umfpack_dl_save_numeric(void* Numeric, char* filename);

    int umfpack_zi_save_numeric(void* Numeric, char* filename);

    UF_long umfpack_zl_save_numeric(void* Numeric, char* filename);

    /* ========================================================================== */
    /* === umfpack_load_numeric ================================================= */
    /* ========================================================================== */

    int umfpack_di_load_numeric(void** Numeric, char* filename);

    UF_long umfpack_dl_load_numeric(void** Numeric, char* filename);

    int umfpack_zi_load_numeric(void** Numeric, char* filename);

    UF_long umfpack_zl_load_numeric(void** Numeric, char* filename);

    /* ========================================================================== */
    /* === umfpack_save_symbolic================================================= */
    /* ========================================================================== */

    int umfpack_di_save_symbolic(void* Symbolic, char* filename);

    UF_long umfpack_dl_save_symbolic(void* Symbolic, char* filename);

    int umfpack_zi_save_symbolic(void* Symbolic, char* filename);

    UF_long umfpack_zl_save_symbolic(void* Symbolic, char* filename);

    /* ========================================================================== */
    /* === umfpack_load_symbolic ================================================ */
    /* ========================================================================== */

    int umfpack_di_load_symbolic(void** Symbolic, char* filename);

    UF_long umfpack_dl_load_symbolic(void** Symbolic, char* filename);

    int umfpack_zi_load_symbolic(void** Symbolic, char* filename);

    UF_long umfpack_zl_load_symbolic(void** Symbolic, char* filename);
    /* ========================================================================== */
    /* === UMFPACK_get_determinant ============================================== */
    /* ========================================================================== */

    /* -------------------------------------------------------------------------- */
    /* UMFPACK Copyright (c) Timothy A. Davis, CISE,                              */
    /* Univ. of Florida.  All Rights Reserved.  See ../Doc/License for License.   */
    /* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
    /* UMFPACK_get_determinant contributed by David Bateman, Motorola, Paris. */
    /* -------------------------------------------------------------------------- */

    int umfpack_di_get_determinant(double* Mx, double* Ex, void* NumericHandle, double* User_Info /*[UMFPACK_INFO]*/
    );

    UF_long umfpack_dl_get_determinant(double* Mx, double* Ex,
            void* NumericHandle, double* User_Info /*[UMFPACK_INFO]*/
    );

    int umfpack_zi_get_determinant(double* Mx, double* Mz, double* Ex,
            void* NumericHandle, double* User_Info /*[UMFPACK_INFO]*/
    );

    UF_long umfpack_zl_get_determinant(double* Mx, double* Mz, double* Ex,
            void* NumericHandle, double* User_Info /*[UMFPACK_INFO]*/
    );

    /* ========================================================================== */
    /* === umfpack_report_status ================================================ */
    /* ========================================================================== */

    void umfpack_di_report_status(/*const*/
            double* Control /*[UMFPACK_CONTROL]*/ , int status);

    void umfpack_dl_report_status(/*const*/
            double* Control /*[UMFPACK_CONTROL]*/ , UF_long status);

    void umfpack_zi_report_status(/*const*/
            double* Control /*[UMFPACK_CONTROL]*/ , int status);

    void umfpack_zl_report_status(/*const*/
            double* Control /*[UMFPACK_CONTROL]*/ , UF_long status);

    /* ========================================================================== */
    /* === umfpack_report_info ================================================== */
    /* ========================================================================== */

    void umfpack_di_report_info(/*const*/
            double* Control /*[UMFPACK_CONTROL]*/ ,/*const*/
            double* Info /*[UMFPACK_INFO]*/
    );

    void umfpack_dl_report_info(/*const*/
            double* Control /*[UMFPACK_CONTROL]*/ ,/*const*/
            double* Info /*[UMFPACK_INFO]*/
    );

    void umfpack_zi_report_info(/*const*/
            double* Control /*[UMFPACK_CONTROL]*/ ,/*const*/
            double* Info /*[UMFPACK_INFO]*/
    );

    void umfpack_zl_report_info(/*const*/
            double* Control /*[UMFPACK_CONTROL]*/ ,/*const*/
            double* Info /*[UMFPACK_INFO]*/
    );

    /* ========================================================================== */
    /* === umfpack_report_control =============================================== */
    /* ========================================================================== */

    void umfpack_di_report_control(/*const*/
            double* Control /*[UMFPACK_CONTROL]*/
    );

    void umfpack_dl_report_control(/*const*/
            double* Control /*[UMFPACK_CONTROL]*/
    );

    void umfpack_zi_report_control(/*const*/
            double* Control /*[UMFPACK_CONTROL]*/
    );

    void umfpack_zl_report_control(/*const*/
            double* Control /*[UMFPACK_CONTROL]*/
    );

    /* ========================================================================== */
    /* === umfpack_report_matrix ================================================ */
    /* ========================================================================== */

    int umfpack_di_report_matrix(int n_row, int n_col,/*const*/
            int* Ap,/*const*/
            int* Ai,/*const*/
            double* Ax,
            int col_form,/*const*/
            double* Control /*[UMFPACK_CONTROL]*/
    );

    UF_long umfpack_dl_report_matrix(UF_long n_row, UF_long n_col,/*const*/
            UF_long* Ap,
            /*const*/
            UF_long* Ai,/*const*/
            double* Ax, UF_long col_form,/*const*/
            double* Control /*[UMFPACK_CONTROL]*/
    );

    int umfpack_zi_report_matrix(int n_row, int n_col,/*const*/
            int* Ap,/*const*/
            int* Ai,/*const*/
            double* Ax, /*const*/ double* Az, int col_form,/*const*/
            double* Control /*[UMFPACK_CONTROL]*/
    );

    UF_long umfpack_zl_report_matrix(UF_long n_row, UF_long n_col,/*const*/
            UF_long* Ap,
            /*const*/
            UF_long* Ai,/*const*/
            double* Ax, /*const*/ double* Az, UF_long col_form,/*const*/
            double* Control /*[UMFPACK_CONTROL]*/
    );

    /* ========================================================================== */
    /* === umfpack_report_triplet =============================================== */
    /* ========================================================================== */

    int umfpack_di_report_triplet(int n_row, int n_col, int nz,/*const*/
            int* Ti,/*const*/
            int* Tj,
            /*const*/
            double* Tx,/*const*/
            double* Control /*[UMFPACK_CONTROL]*/
    );

    UF_long umfpack_dl_report_triplet(UF_long n_row, UF_long n_col, UF_long nz,
            /*const*/
            UF_long* Ti,/*const*/
            UF_long* Tj,/*const*/
            double* Tx,/*const*/
            double* Control /*[UMFPACK_CONTROL]*/
    );

    int umfpack_zi_report_triplet(int n_row, int n_col, int nz,/*const*/
            int* Ti,/*const*/
            int* Tj,
            /*const*/
            double* Tx, /*const*/ double* Tz,/*const*/
            double* Control /*[UMFPACK_CONTROL]*/
    );

    UF_long umfpack_zl_report_triplet(UF_long n_row, UF_long n_col, UF_long nz,
            /*const*/
            UF_long* Ti,/*const*/
            UF_long* Tj,/*const*/
            double* Tx, /*const*/ double* Tz,/*const*/
            double* Control /*[UMFPACK_CONTROL]*/
    );

    /* ========================================================================== */
    /* === umfpack_report_vector ================================================ */
    /* ========================================================================== */

    int umfpack_di_report_vector(int n,/*const*/
            double* X,/*const*/
            double* Control /*[UMFPACK_CONTROL]*/
    );

    UF_long umfpack_dl_report_vector(UF_long n,/*const*/
            double* X,/*const*/
            double* Control /*[UMFPACK_CONTROL]*/
    );

    int umfpack_zi_report_vector(int n,/*const*/
            double* Xx, /*const*/ double* Xz,/*const*/
            double* Control /*[UMFPACK_CONTROL]*/
    );

    UF_long umfpack_zl_report_vector(UF_long n,/*const*/
            double* Xx, /*const*/ double* Xz,/*const*/
            double* Control /*[UMFPACK_CONTROL]*/
    );

    /* ========================================================================== */
    /* === umfpack_report_symbolic ============================================== */
    /* ========================================================================== */

    int umfpack_di_report_symbolic(void* Symbolic,/*const*/
            double* Control /*[UMFPACK_CONTROL]*/
    );

    UF_long umfpack_dl_report_symbolic(void* Symbolic,/*const*/
            double* Control /*[UMFPACK_CONTROL]*/
    );

    int umfpack_zi_report_symbolic(void* Symbolic,/*const*/
            double* Control /*[UMFPACK_CONTROL]*/
    );

    UF_long umfpack_zl_report_symbolic(void* Symbolic,/*const*/
            double* Control /*[UMFPACK_CONTROL]*/
    );

    /* ========================================================================== */
    /* === umfpack_report_numeric =============================================== */
    /* ========================================================================== */

    int umfpack_di_report_numeric(void* Numeric,/*const*/
            double* Control /*[UMFPACK_CONTROL]*/
    );

    UF_long umfpack_dl_report_numeric(void* Numeric,/*const*/
            double* Control /*[UMFPACK_CONTROL]*/
    );

    int umfpack_zi_report_numeric(void* Numeric,/*const*/
            double* Control /*[UMFPACK_CONTROL]*/
    );

    UF_long umfpack_zl_report_numeric(void* Numeric,/*const*/
            double* Control /*[UMFPACK_CONTROL]*/
    );

    /* ========================================================================== */
    /* === umfpack_report_perm ================================================== */
    /* ========================================================================== */

    int umfpack_di_report_perm(int np,/*const*/
            int* Perm,/*const*/
            double* Control /*[UMFPACK_CONTROL]*/
    );

    UF_long umfpack_dl_report_perm(UF_long np,/*const*/
            UF_long* Perm,/*const*/
            double* Control /*[UMFPACK_CONTROL]*/
    );

    int umfpack_zi_report_perm(int np,/*const*/
            int* Perm,/*const*/
            double* Control /*[UMFPACK_CONTROL]*/
    );

    UF_long umfpack_zl_report_perm(UF_long np,/*const*/
            UF_long* Perm,/*const*/
            double* Control /*[UMFPACK_CONTROL]*/
    );

    /* ========================================================================== */
    /* === umfpack_timer ======================================================== */
    /* ========================================================================== */

    double umfpack_timer();

    /* ========================================================================== */
    /* === umfpack_tictoc ======================================================= */
    /* ========================================================================== */

    void umfpack_tic(double* stats /*[2]*/ );

    void umfpack_toc(double* stats /*[2]*/ );

} // end extern(C)

/* UMFPACK Version 4.5 and later will include the following definitions.
 * As an example, to test if the version you are using is 4.5 or later:
 *
 * #ifdef UMFPACK_VER
 *	if (UMFPACK_VER >= UMFPACK_VER_CODE (4,5)) ...
 * #endif
 *
 * This also works during compile-time:
 *
 *	#if defined(UMFPACK_VER) && (UMFPACK >= UMFPACK_VER_CODE (4,5))
 *	    printf ("This is version 4.5 or later\n") ;
 *	#else
 *	    printf ("This is an early version\n") ;
 *	#endif
 *
 * Versions 4.4 and earlier of UMFPACK do not include a #define'd version
 * number, although they do include the UMFPACK_VERSION string, defined
 * above.
 */

const char[] UMFPACK_VERSION = "UMFPACK V5.1.0 (May 31, 2007)";
const char[] UMFPACK_DATE = "May 31, 2007";

/*
UMFPACK: Copyright 1995-2006 by Timothy A. Davis.  All Rights Reserved.
UMFPACK is available under alternate licenses, contact T. Davis for details.

UMFPACK License:

    Your use or distribution of UMFPACK or any modified version of
    UMFPACK implies that you agree to this License.

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with this library; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301
    USA

    Permission is hereby granted to use or copy this program under the
    terms of the GNU LGPL, provided that the Copyright, this License,
    and the Availability of the original version is retained on all copies.
    User documentation of any code that uses this code or any modified
    version of this code must cite the Copyright, this License, the
    Availability note, and "Used by permission." Permission to modify
    the code and to distribute modified code is granted, provided the
    Copyright, this License, and the Availability note are retained,
    and a notice that the code was modified is included.

Availability:

    http://www.cise.ufl.edu/research/sparse/umfpack
*/

enum
{
    UMFPACK_MAIN_VERSION = 5,
    UMFPACK_SUB_VERSION = 1,
    UMFPACK_SUBSUB_VERSION = 0,
    UMFPACK_VER = 5001,
}

/* -------------------------------------------------------------------------- */
/* contents of Info */
/* -------------------------------------------------------------------------- */

/* Note that umfpack_report.m must coincide with these definitions.  S is
 * the submatrix of A after removing row/col singletons and empty rows/cols. */

enum
{
    /* returned by all routines that use Info: */
    UMFPACK_STATUS = 0, /* UMFPACK_OK, or other result */
    UMFPACK_NROW = 1, /* n_row input value */
    UMFPACK_NCOL = 16, /* n_col input value */
    UMFPACK_NZ = 2, /* # of entries in A */

    /* computed in UMFPACK_*symbolic and UMFPACK_numeric: */
    UMFPACK_SIZE_OF_UNIT = 3, /* sizeof (Unit) */

    /* computed in UMFPACK_*symbolic: */
    UMFPACK_SIZE_OF_INT = 4, /* sizeof (int) */
    UMFPACK_SIZE_OF_LONG = 5, /* sizeof (UF_long) */
    UMFPACK_SIZE_OF_POINTER = 6, /* sizeof (void *) */
    UMFPACK_SIZE_OF_ENTRY = 7, /* sizeof (Entry), real or complex */
    UMFPACK_NDENSE_ROW = 8, /* number of dense rows */
    UMFPACK_NEMPTY_ROW = 9, /* number of empty rows */
    UMFPACK_NDENSE_COL = 10, /* number of dense rows */
    UMFPACK_NEMPTY_COL = 11, /* number of empty rows */
    UMFPACK_SYMBOLIC_DEFRAG = 12, /* # of memory compactions */
    UMFPACK_SYMBOLIC_PEAK_MEMORY = 13, /* memory used by symbolic analysis */
    UMFPACK_SYMBOLIC_SIZE = 14, /* size of Symbolic object, in Units */
    UMFPACK_SYMBOLIC_TIME = 15, /* time (sec.) for symbolic analysis */
    UMFPACK_SYMBOLIC_WALLTIME = 17, /* wall clock time for sym. analysis */
    UMFPACK_STRATEGY_USED = 18, /* strategy used: sym, unsym, 2by2 */
    UMFPACK_ORDERING_USED = 19, /* ordering used: colamd, amd, given */
    UMFPACK_QFIXED = 31, /* whether Q is fixed or refined */
    UMFPACK_DIAG_PREFERRED = 32, /* whether diagonal pivoting attempted*/
    UMFPACK_PATTERN_SYMMETRY = 33, /* symmetry of pattern of S */
    UMFPACK_NZ_A_PLUS_AT = 34, /* nnz (S+S'), excl. diagonal */
    UMFPACK_NZDIAG = 35, /* nnz (diag (S)) */

    /* AMD statistics, computed in UMFPACK_*symbolic: */
    UMFPACK_SYMMETRIC_LUNZ = 36, /* nz in L+U, if AMD ordering used */
    UMFPACK_SYMMETRIC_FLOPS = 37, /* flops for LU, if AMD ordering used */
    UMFPACK_SYMMETRIC_NDENSE = 38, /* # of "dense" rows/cols in S+S' */
    UMFPACK_SYMMETRIC_DMAX = 39, /* max nz in cols of L, for AMD */

    /* statistics for 2-by-2 strategy */
    UMFPACK_2BY2_NWEAK = 51, /* number of weak diagonal entries*/
    UMFPACK_2BY2_UNMATCHED = 52, /* # of weak diagonals not matched*/
    UMFPACK_2BY2_PATTERN_SYMMETRY = 53, /* symmetry of pattern of P*S */
    UMFPACK_2BY2_NZ_PA_PLUS_PAT = 54, /* nz in PS+(PS)' */
    UMFPACK_2BY2_NZDIAG = 55, /* nz on diagonal of PS+(PS)' */

    /* statistcs for singleton pruning */
    UMFPACK_COL_SINGLETONS = 56, /* # of column singletons */
    UMFPACK_ROW_SINGLETONS = 57, /* # of row singletons */
    UMFPACK_N2 = 58, /* size of S */
    UMFPACK_S_SYMMETRIC = 59, /* 1 if S square and symmetricly perm.*/

    /* estimates computed in UMFPACK_*symbolic: */
    UMFPACK_NUMERIC_SIZE_ESTIMATE = 20, /* final size of Numeric->Memory */
    UMFPACK_PEAK_MEMORY_ESTIMATE = 21, /* for symbolic & numeric */
    UMFPACK_FLOPS_ESTIMATE = 22, /* flop count */
    UMFPACK_LNZ_ESTIMATE = 23, /* nz in L, incl. diagonal */
    UMFPACK_UNZ_ESTIMATE = 24, /* nz in U, incl. diagonal */
    UMFPACK_VARIABLE_INIT_ESTIMATE = 25, /* initial size of Numeric->Memory*/
    UMFPACK_VARIABLE_PEAK_ESTIMATE = 26, /* peak size of Numeric->Memory */
    UMFPACK_VARIABLE_FINAL_ESTIMATE = 27, /* final size of Numeric->Memory */
    UMFPACK_MAX_FRONT_SIZE_ESTIMATE = 28, /* max frontal matrix size */
    UMFPACK_MAX_FRONT_NROWS_ESTIMATE = 29, /* max # rows in any front */
    UMFPACK_MAX_FRONT_NCOLS_ESTIMATE = 30, /* max # columns in any front */

    /* exact values, (estimates shown above) computed in UMFPACK_numeric: */
    UMFPACK_NUMERIC_SIZE = 40, /* final size of Numeric->Memory */
    UMFPACK_PEAK_MEMORY = 41, /* for symbolic & numeric */
    UMFPACK_FLOPS = 42, /* flop count */
    UMFPACK_LNZ = 43, /* nz in L, incl. diagonal */
    UMFPACK_UNZ = 44, /* nz in U, incl. diagonal */
    UMFPACK_VARIABLE_INIT = 45, /* initial size of Numeric->Memory*/
    UMFPACK_VARIABLE_PEAK = 46, /* peak size of Numeric->Memory */
    UMFPACK_VARIABLE_FINAL = 47, /* final size of Numeric->Memory */
    UMFPACK_MAX_FRONT_SIZE = 48, /* max frontal matrix size */
    UMFPACK_MAX_FRONT_NROWS = 49, /* max # rows in any front */
    UMFPACK_MAX_FRONT_NCOLS = 50, /* max # columns in any front */

    /* computed in UMFPACK_numeric: */
    UMFPACK_NUMERIC_DEFRAG = 60, /* # of garbage collections */
    UMFPACK_NUMERIC_REALLOC = 61, /* # of memory reallocations */
    UMFPACK_NUMERIC_COSTLY_REALLOC = 62, /* # of costlly memory realloc's */
    UMFPACK_COMPRESSED_PATTERN = 63, /* # of integers in LU pattern */
    UMFPACK_LU_ENTRIES = 64, /* # of reals in LU factors */
    UMFPACK_NUMERIC_TIME = 65, /* numeric factorization time */
    UMFPACK_UDIAG_NZ = 66, /* nz on diagonal of U */
    UMFPACK_RCOND = 67, /* est. reciprocal condition # */
    UMFPACK_WAS_SCALED = 68, /* none, max row, or sum row */
    UMFPACK_RSMIN = 69, /* min (max row) or min (sum row) */
    UMFPACK_RSMAX = 70, /* max (max row) or max (sum row) */
    UMFPACK_UMIN = 71, /* min abs diagonal entry of U */
    UMFPACK_UMAX = 72, /* max abs diagonal entry of U */
    UMFPACK_ALLOC_INIT_USED = 73, /* alloc_init parameter used */
    UMFPACK_FORCED_UPDATES = 74, /* # of forced updates */
    UMFPACK_NUMERIC_WALLTIME = 75, /* numeric wall clock time */
    UMFPACK_NOFF_DIAG = 76, /* number of off-diagonal pivots */

    UMFPACK_ALL_LNZ = 77, /* nz in L, if no dropped entries */
    UMFPACK_ALL_UNZ = 78, /* nz in U, if no dropped entries */
    UMFPACK_NZDROPPED = 79, /* # of dropped small entries */

    /* computed in UMFPACK_solve: */
    UMFPACK_IR_TAKEN = 80, /* # of iterative refinement steps taken */
    UMFPACK_IR_ATTEMPTED = 81, /* # of iter. refinement steps attempted */
    UMFPACK_OMEGA1 = 82, /* omega1, sparse backward error estimate */
    UMFPACK_OMEGA2 = 83, /* omega2, sparse backward error estimate */
    UMFPACK_SOLVE_FLOPS = 84, /* flop count for solve */
    UMFPACK_SOLVE_TIME = 85, /* solve time (seconds) */
    UMFPACK_SOLVE_WALLTIME = 86, /* solve time (wall clock, seconds) */

    /* Info [87, 88, 89] unused */

}
/* Unused parts of Info may be used in future versions of UMFPACK. */

/* -------------------------------------------------------------------------- */

enum
{
    /* Info [UMFPACK_ORDERING_USED] is one of the following: */
    UMFPACK_ORDERING_COLAMD = 0, /* COLAMD(A) */
    UMFPACK_ORDERING_AMD = 1, /* AMD(A+A') */
    UMFPACK_ORDERING_GIVEN = 2, /* Q is provided on input */

}

/* -------------------------------------------------------------------------- */
/* contents of Control */
/* -------------------------------------------------------------------------- */

enum
{
    /* used in all UMFPACK_report_* routines: */
    UMFPACK_PRL = 0, /* print level */

    /* used in UMFPACK_*symbolic only: */
    UMFPACK_DENSE_ROW = 1, /* dense row parameter */
    UMFPACK_DENSE_COL = 2, /* dense col parameter */
    UMFPACK_BLOCK_SIZE = 4, /* BLAS-3 block size */
    UMFPACK_STRATEGY = 5, /* auto, symmetric, unsym., or 2by2 */
    UMFPACK_2BY2_TOLERANCE = 12, /* 2-by-2 pivot tolerance */
    UMFPACK_FIXQ = 13, /* -1: no fixQ, 0: default, 1: fixQ */
    UMFPACK_AMD_DENSE = 14, /* for AMD ordering */
    UMFPACK_AGGRESSIVE = 19, /* whether or not to use aggressive
					 * absorption in AMD and COLAMD */

    /* used in UMFPACK_numeric only: */
    UMFPACK_PIVOT_TOLERANCE = 3, /* threshold partial pivoting setting */
    UMFPACK_ALLOC_INIT = 6, /* initial allocation ratio */
    UMFPACK_SYM_PIVOT_TOLERANCE = 15, /* threshold, only for diag. entries */
    UMFPACK_SCALE = 16, /* what row scaling to do */
    UMFPACK_FRONT_ALLOC_INIT = 17, /* frontal matrix allocation ratio */
    UMFPACK_DROPTOL = 18, /* drop tolerance for entries in L,U */

    /* used in UMFPACK_*solve only: */
    UMFPACK_IRSTEP = 7, /* max # of iterative refinements */

    /* compile-time settings - Control [8..11] cannot be changed at run time: */
    UMFPACK_COMPILED_WITH_BLAS = 8, /* uses the BLAS */
    UMFPACK_COMPILED_FOR_MATLAB = 9, /* 1 if MATLAB mexFunction, etc. */
    UMFPACK_COMPILED_WITH_GETRUSAGE = 10, /* uses getrusage timer, or not */
    UMFPACK_COMPILED_IN_DEBUG_MODE = 11, /* debugging enabled (very slow!) */

    /* -------------------------------------------------------------------------- */

    /* Control [UMFPACK_STRATEGY] is one of the following: */
    UMFPACK_STRATEGY_AUTO = 0, /* use sym. or unsym. strategy */
    UMFPACK_STRATEGY_UNSYMMETRIC = 1, /* COLAMD(A), coletree postorder,
                                       not prefer diag*/
    UMFPACK_STRATEGY_2BY2 = 2, /* AMD(PA+PA'), no coletree postorder,
                                   prefer diag(PA) where P is pseudo
                                   max transversal */
    UMFPACK_STRATEGY_SYMMETRIC = 3, /* AMD(A+A'), no coletree postorder,
                                   prefer diagonal */

    /* Control [UMFPACK_SCALE] is one of the following: */
    UMFPACK_SCALE_NONE = 0, /* no scaling */
    UMFPACK_SCALE_SUM = 1, /* default: divide each row by sum (abs (row))*/
    UMFPACK_SCALE_MAX = 2, /* divide each row by max (abs (row)) */



}
/* -------------------------------------------------------------------------- */
/* default values of Control: */
/* -------------------------------------------------------------------------- */

enum
{
    UMFPACK_DEFAULT_PRL = 1,
    UMFPACK_DEFAULT_BLOCK_SIZE = 32,
    UMFPACK_DEFAULT_IRSTEP = 2,
    UMFPACK_DEFAULT_SCALE = UMFPACK_SCALE_SUM,
    UMFPACK_DEFAULT_STRATEGY = UMFPACK_STRATEGY_AUTO,
    UMFPACK_DEFAULT_FIXQ = 0,
    UMFPACK_DEFAULT_AGGRESSIVE = 1,
}

const double UMFPACK_DEFAULT_AMD_DENSE = 10.0; //=AMD_DEFAULT_DENSE, 
const double UMFPACK_DEFAULT_DENSE_ROW = 0.2;
const double UMFPACK_DEFAULT_DENSE_COL = 0.2;
const double UMFPACK_DEFAULT_PIVOT_TOLERANCE = 0.1;
const double UMFPACK_DEFAULT_2BY2_TOLERANCE = 0.01;
const double UMFPACK_DEFAULT_SYM_PIVOT_TOLERANCE = 0.001;
const double UMFPACK_DEFAULT_ALLOC_INIT = 0.7;
const double UMFPACK_DEFAULT_FRONT_ALLOC_INIT = 0.5;
const double UMFPACK_DEFAULT_DROPTOL = 0;

/* default values of Control may change in future versions of UMFPACK. */

/* -------------------------------------------------------------------------- */
/* status codes */
/* -------------------------------------------------------------------------- */

enum : int
{
    UMFPACK_OK = (0),

    /* status > 0 means a warning, but the method was successful anyway. */
    /* A Symbolic or Numeric object was still created. */
    UMFPACK_WARNING_singular_matrix = (1),

    /* The following warnings were added in umfpack_*_get_determinant */
    UMFPACK_WARNING_determinant_underflow = (2),
    UMFPACK_WARNING_determinant_overflow = (3),

    /* status < 0 means an error, and the method was not successful. */
    /* No Symbolic of Numeric object was created. */
    UMFPACK_ERROR_out_of_memory = (-1),
    UMFPACK_ERROR_invalid_Numeric_object = (-3),
    UMFPACK_ERROR_invalid_Symbolic_object = (-4),
    UMFPACK_ERROR_argument_missing = (-5),
    UMFPACK_ERROR_n_nonpositive = (-6),
    UMFPACK_ERROR_invalid_matrix = (-8),
    UMFPACK_ERROR_different_pattern = (-11),
    UMFPACK_ERROR_invalid_system = (-13),
    UMFPACK_ERROR_invalid_permutation = (-15),
    UMFPACK_ERROR_internal_error = (-911), /* yes, call me if you get this! */
    UMFPACK_ERROR_file_IO = (-17),
}
/* -------------------------------------------------------------------------- */
/* solve codes */
/* -------------------------------------------------------------------------- */

/* Solve the system ( )x=b, where ( ) is defined below.  "t" refers to the */
/* linear algebraic transpose (complex conjugate if A is complex), or the (') */
/* operator in MATLAB.  "at" refers to the array transpose, or the (.') */
/* operator in MATLAB. */
enum
{
    UMFPACK_A = (0), /* Ax=b    */
    UMFPACK_At = (1), /* A'x=b   */
    UMFPACK_Aat = (2), /* A.'x=b  */

    UMFPACK_Pt_L = (3), /* P'Lx=b  */
    UMFPACK_L = (4), /* Lx=b    */
    UMFPACK_Lt_P = (5), /* L'Px=b  */
    UMFPACK_Lat_P = (6), /* L.'Px=b */
    UMFPACK_Lt = (7), /* L'x=b   */
    UMFPACK_Lat = (8), /* L.'x=b  */

    UMFPACK_U_Qt = (9), /* UQ'x=b  */
    UMFPACK_U = (10), /* Ux=b    */
    UMFPACK_Q_Ut = (11), /* QU'x=b  */
    UMFPACK_Q_Uat = (12), /* QU.'x=b */
    UMFPACK_Ut = (13), /* U'x=b   */
    UMFPACK_Uat = (14), /* U.'x=b  */



} // end enum

/* -------------------------------------------------------------------------- */

/* Integer constants are used for status and solve codes instead of enum */
/* to make it easier for a Fortran code to call UMFPACK. */

//--- Emacs setup ---
// Local Variables:
// c-basic-offset: 4
// indent-tabs-mode: nil
// End:
