#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <array>
#include <algorithm>
#include <iterator>
#include <math.h>
#include <string.h>
#include <assert.h>

#include <omp.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_fermi_dirac.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_result.h>
#include <gsl/gsl_integration.h>

#define H5_USE_16_API
#include "hdf5.h"

#include "eos.h"
#include "constants.h"
#include "nuclear.h"
#include "physics.h"
#include "abundances.h"

#define max(x, y) ((x)>(y)?(x):(y))
#define min(x, y) ((x)<(y)?(x):(y))
