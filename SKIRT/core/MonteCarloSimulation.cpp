/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "MonteCarloSimulation.hpp"
#include "FatalError.hpp"
#include "Log.hpp"
#include "Random.hpp"
#include "StringUtils.hpp"
#include "System.hpp"
#include "TimeLogger.hpp"

// included for testing purposes
#include "FilePaths.hpp"

////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::setupSelfBefore()
{
    Simulation::setupSelfBefore();
}

////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::setEmulationMode()
{
    _emulationMode = true;
    _numPackages = 0;
}

////////////////////////////////////////////////////////////////////

bool MonteCarloSimulation::emulationMode()
{
    return _emulationMode;
}

////////////////////////////////////////////////////////////////////

int MonteCarloSimulation::dimension() const
{
    return 0;
}

////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::runSelf()
{
    TimeLogger logger(log(), "the test phase");

    for (int k = 1; k<=9; ++k)
    {
        TimeLogger logger(log(), "big table " + std::to_string(k));

        auto map = System::acquireMemoryMap(FilePaths::resource("BigTable"+std::to_string(k)+".stab"));
        if (!map.first) throw FATALERROR("Cannot acquire memory map");
        auto start = static_cast<const char*>(map.first);
        auto end = static_cast<const char*>(map.first) + map.second;

        log()->info(string(start, 7) + " ... " + string(end-8, 7));
        size_t length = *(static_cast<const size_t*>(map.first) + 6);
        log()->info("length: " + std::to_string(length));

        const double* x = static_cast<const double*>(map.first) + 7;
        const double* y = static_cast<const double*>(map.first) + 7 + length + 4;
        log()->info("x: " + StringUtils::toString(x[1]) + "    y: " + StringUtils::toString(y[1]));

        double xsum = 0;
        double ysum = 0;
        for (size_t i = 0; i<length; ++i)
        {
            xsum += x[i];
            ysum += y[i];
        }
        log()->warning("xsum: " + StringUtils::toString(xsum) + "    ysum: " + StringUtils::toString(ysum));
    }
}

////////////////////////////////////////////////////////////////////
