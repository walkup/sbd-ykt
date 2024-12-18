#ifndef SBD_HCBOSON_OUT_OF_PLACE_FUNC_INC_ALL_H
#define SBD_HCBOSON_OUT_OF_PLACE_FUNC_INC_ALL_H

#include "sbd/hcboson/out_of_place_func/arithmetic.h"
#ifdef SBD_USE_NONBLOCKING
#include "sbd/hcboson/out_of_place_func/mult_nb.h"
#else
#include "sbd/hcboson/out_of_place_func/mult_bc.h"
#endif
#include "sbd/hcboson/out_of_place_func/lanczos.h"

#endif
