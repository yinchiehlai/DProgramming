/**
 * Simple binding for multi-thread FFTW3  (http://www.fftw.org/) .
 */

module addon.fftw_thread;
import addon.tool;

pragma(lib, "fftw3_threads");
pragma(lib, "fftw3");

extern (C)
{
    void* fftw_plan_dft_1d(int n, Cdouble* vin, Cdouble* vout, int FFTW_FORWARD, int ffttype);
    void* fftw_plan_dft_2d(int dim0, int dim1, Cdouble* vin, Cdouble* vout,
            int FFTW_FORWARD, int ffttype);
    void fftw_execute(void* p1);
    void fftw_destroy_plan(void* p1);
    int fftw_init_threads();
    void fftw_plan_with_nthreads(int n);
    void fftw_cleanup_threads();
}

const size_t FFTW_MEASURE = (0U);
const size_t FFTW_ESTIMATE = (1U << 6);
const int FFTW_FORWARD = -1;
const int FFTW_BACKWARD = +1;

void* fftw_plan(Cdouble[] vin, Cdouble[] vout)
{
    return fftw_plan_dft_1d(cast(int) vin.length, vin.ptr, vout.ptr, FFTW_FORWARD, FFTW_ESTIMATE);
}


void* ifftw_plan(Cdouble[] vin, Cdouble[] vout)
{
    return fftw_plan_dft_1d(cast(int) vin.length, vin.ptr, vout.ptr, FFTW_BACKWARD, FFTW_ESTIMATE);
}

void* fftw_2d_plan(CIdx vin, CIdx vout)
{
    return fftw_plan_dft_2d(cast(int) vin.dim[0], cast(int) vin.dim[1], vin.srg.ptr, vout.srg.ptr, FFTW_FORWARD, FFTW_ESTIMATE);
}


void* ifftw_2d_plan(CIdx vin, CIdx vout)
{
    return fftw_plan_dft_2d(cast(int) vin.dim[0], cast(int) vin.dim[1], vin.srg.ptr, vout.srg.ptr, FFTW_BACKWARD, FFTW_ESTIMATE);
}
