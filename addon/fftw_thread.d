/**
 * Simple binding for multi-thread FFTW3  (http://www.fftw.org/) .
 */

module addon.fftw_thread;
import addon.tool;

pragma(lib, "fftw3_threads");

extern (C)
{
    void* fftw_plan_dft_1d(int n, Cdouble* vin, Cdouble* vout, int FFTW_FORWARD, int ffttype);
    void* fftw_plan_dft_2d(int dim0, int dim1, Cdouble* vin, Cdouble* vout,
            int FFTW_FORWARD, int ffttype);
    void fftw_execute(void* p1);
    void fftw_destroy_plan(void* p1);
    void fftw_init_threads();
    void fftw_plan_with_nthreads(int n);
    void fftw_cleanup_threads();
}

const size_t FFTW_MEASURE = (0U);
const size_t FFTW_ESTIMATE = (1U << 6);
const int FFTW_FORWARD = -1;
const int FFTW_BACKWARD = +1;

void fft(Cdouble[] vin, Cdouble[] vout)
{
    void* p1 = fftw_plan_dft_1d(cast(int) vin.length, vin.ptr, vout.ptr,
            FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p1);
    fftw_end(p1);
}

void* fftw_plan(Cdouble[] vin, Cdouble[] vout)
{
    return fftw_plan_dft_1d(cast(int) vin.length, vin.ptr, vout.ptr, FFTW_FORWARD, FFTW_ESTIMATE);
}

void ifft(Cdouble[] vin, Cdouble[] vout)
{
    void* p1 = fftw_plan_dft_1d(cast(int) vin.length, vin.ptr, vout.ptr,
            FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(p1);
    fftw_end(p1);
}

void* ifftw_plan(Cdouble[] vin, Cdouble[] vout)
{
    return fftw_plan_dft_1d(cast(int) vin.length, vin.ptr, vout.ptr, FFTW_BACKWARD, FFTW_ESTIMATE);
}

void fftw_end(void* p1)
{
    fftw_destroy_plan(p1);
}
