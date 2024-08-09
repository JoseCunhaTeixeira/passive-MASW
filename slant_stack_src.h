#pragma once
#ifndef SLANT_STACK_SRC_H
#define SLANT_STACK_SRC_H

#include <vector>
#include <complex>

std::vector<std::vector<double> > makeFV_src(const std::vector<std::vector<double> >& XT, double si, const std::vector<double>& offset, double vmin, double vmax, double dv, double fmax);

#endif 
