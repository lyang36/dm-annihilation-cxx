%module vl2read

%{
#include <string>
#include <cstring>
#include <cmath>
#include <rpc/types.h>
#include <rpc/xdr.h>
#include <vector>
#include "../tipsydefs.h"
#include "../datareader.h"
#include "vl2read.h"
%}

%include <typemaps.i>
%include <carrays.i>
%include <std_string.i>
%include <std_vector.i>
%include "array_util.i"


%include "vl2read.h"

#double array, use this one to wap the array
%array_class(double, dArray);
%array_class(float, fArray);
