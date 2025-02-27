/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef FSPSSED_HPP
#define FSPSSED_HPP

#include "FamilySED.hpp"

////////////////////////////////////////////////////////////////////

/** FSPSSED is a class that represents spectral energy distributions of simple stellar populations
    (SSPs), parameterized on metallicity and age, generated by the FSPS code presented by Conroy,
    Gunn, & White (2009, ApJ, 699, 486) and Conroy & Gunn (2010, ApJ, 712, 833), and assuming a
    Chabrier, Kroupa or Salpeter initial mass function. See the FSPSSEDFamily class for more
    information. */
class FSPSSED : public FamilySED
{
    /** The enumeration type indicating the assumed initial mass function (IMF). */
    ENUM_DEF(IMF, Chabrier, Kroupa, Salpeter)
    ENUM_VAL(IMF, Chabrier, "Chabrier IMF")
    ENUM_VAL(IMF, Kroupa, "Kroupa IMF")
    ENUM_VAL(IMF, Salpeter, "Salpeter IMF")
    ENUM_END()

    ITEM_CONCRETE(FSPSSED, FamilySED, "an FSPS simple stellar population SED")

    PROPERTY_ENUM(imf, IMF, "the assumed initial mass function")
        ATTRIBUTE_DEFAULT_VALUE(imf, "Chabrier")

    PROPERTY_DOUBLE(metallicity, "the metallicity of the SSP")
        ATTRIBUTE_MIN_VALUE(metallicity, "[0.0001")
        ATTRIBUTE_MAX_VALUE(metallicity, "0.05]")
        ATTRIBUTE_DEFAULT_VALUE(metallicity, "0.02")

    PROPERTY_DOUBLE(age, "the age of the SSP")
        ATTRIBUTE_QUANTITY(age, "time")
        ATTRIBUTE_MIN_VALUE(age, "[0 Gyr")
        ATTRIBUTE_MAX_VALUE(age, "20 Gyr]")
        ATTRIBUTE_DEFAULT_VALUE(age, "5 Gyr")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function returns a newly created SEDFamily object (which is already hooked into the
        simulation item hierachy so it will be automatically deleted) and stores the parameters for
        the specific %SED configured by the user in the specified array. */
    const SEDFamily* getFamilyAndParameters(Array& parameters) override;
};

////////////////////////////////////////////////////////////////////

#endif
