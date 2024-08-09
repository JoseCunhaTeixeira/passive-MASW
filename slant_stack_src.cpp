
#include "slant_stack_src.h"
#include <iostream>
#include <cmath>
#include <fftw3.h>
#include <cstdio>
#include <omp.h>

#define NB_THREADS 4

/**
 * @author J. Cunha Teixeira
 * @author B. Decker
 * 
 * Constructs a FV dispersion diagram
 * @param XT data
 * @param si sampling interval in seconds
 * @param offsets offsets in meter
 * @param vmin minimal velocity to scan in m/s
 * @param vmax maximal velocity to scan in m/s
 * @param dv velocity step in m/s
 * @param fmax maximum frequency computed
 * @return a matrix FV of size vs.size() * fs.size() representing a dispersing plot
*/
// Cunha Teixeira - Implemented in Python
// JUL2024 - Decker - Traduced in C++
std::vector<std::vector<double> > makeFV_src(const std::vector<std::vector<double> >& XT, 
                                            double si, 
                                            const std::vector<double>& offset, 
                                            double vmin, 
                                            double vmax, 
                                            double dv, 
                                            double fmax) {
    omp_set_num_threads(NB_THREADS);
    int Nt = XT[0].size();
    int Nx = XT.size();

    // Transformation de Fourier avec FFTW
    fftw_complex* in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * Nt);
    fftw_complex* out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * Nt);
    fftw_plan p = fftw_plan_dft_1d(Nt, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    
    std::vector<std::vector<std::complex<double> > > XF(Nx, std::vector<std::complex<double> >(Nt));
    #pragma omp parallel for
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Nt; ++j) {
            in[j][0] = XT[i][j]; // partie réelle
            in[j][1] = 0.0; // partie imaginaire
        }
        fftw_execute(p);
        for (int j = 0; j < Nt; ++j) {
            XF[i][j] = std::complex<double>(out[j][0], out[j][1]);
        }
    }
    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(out);
    // Calcul de l'axe des fréquences
    std::vector<double> fs(Nt);
    for (int i = 0; i < Nt; ++i) {
        fs[i] = i / (si * Nt);
    }
    
    int imax = 0;
    for (int i = 0; i < Nt; ++i) {
        if (fs[i] > fmax) {
            imax = i;
            break;
        }
    }
    fs.resize(imax);
    // Tronquer XF à la plage de fréquences d'intérêt
    for (int i = 0; i < Nx; ++i) {
        XF[i].resize(imax);
    }
    
    // Axe des vitesses
    std::vector<double> vs;
    for (double v = vmin; v <= vmax; v += dv) {
        vs.push_back(v);
    }

    std::vector<std::vector<double> > FV(vs.size(), std::vector<double>(fs.size(), 0.0));

    #pragma omp parallel for
    for (size_t i = 0; i < vs.size(); ++i) {
        double v = vs[i];
        for (size_t j = 0; j < fs.size(); ++j) {
            double f = fs[j];
            std::complex<double> sum_exp(0.0, 0.0);
            for (size_t k = 0; k < offset.size(); ++k) {
                double dphi = 2 * M_PI * offset[k] * f / v;
                std::complex<double> exp_term = std::polar(1., dphi) * (XF[k][j] / abs(XF[k][j]));
                sum_exp += exp_term;
            }
            FV[i][j] = std::abs(sum_exp);
        }
    }
    return FV;
}


