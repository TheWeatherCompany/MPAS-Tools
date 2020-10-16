#ifndef _REMAPPERBASE_H
#define _REMAPPERBASE_H

#include <typeinfo>

class RemapperBase {
public:
	RemapperBase();
	~RemapperBase();

    enum interp_type
    {
        barycentric,
        nearest,
        nearest_land,
        nearest_water,
        nearest_samelandmask
    };

    virtual void remap(const std::type_info& t, int ndims, interp_type interp, void *dst, void *src) = 0;
};
#endif
