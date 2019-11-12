/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "HDF5InFile.hpp"
#include "FatalError.hpp"
#include "FilePaths.hpp"
#include "Log.hpp"
#include "StringUtils.hpp"
#include "System.hpp"
#include "Units.hpp"
#include <exception>
#include <regex>
#include <sstream>



////////////////////////////////////////////////////////////////////

HDF5InFile::HDF5InFile(const SimulationItem* item, string filename, string description)
{
    // open the file
    string filepath = item->find<FilePaths>()->input(filename);
    _in = new HF::File(filepath,HF::File::ReadOnly);

    // figure out the columns that we have
    std::vector<std::string> object_names = _in->listObjectNames();
    size_t index=0;
    for(int i=0;i<object_names.size();i++)
    {
        if(_in->getObjectType(object_names[i])==HF::ObjectType::Dataset)
        {
            std::string name = object_names[i];
            std::cout << name << " is a dataset"  << std::endl;
            HF::DataSet ds = _in->getDataSet(name);
            HF::Attribute att = ds.getAttribute("unit");
            std::string unit;
            att.read(unit);
            _colv.emplace_back();
            _colv.back().physColIndex = ++index;
            _colv.back().unit = unit;	
            _colv.back().title = name;
        }
    }
    // remember the units system and the logger
    _units = item->find<Units>();
    _log = item->find<Log>();

    // log "reading file" message
    _log->info( item->typeAndName() + " reads " + description + " from text file " + filepath + "...");

    _hasFileInfo = !_colv.empty();
}

////////////////////////////////////////////////////////////////////

void HDF5InFile::close()
{
    if(_in)
        delete _in;
}

//////////////////////////////////////////////////////////////////////
//
HDF5InFile::~HDF5InFile()
{
    close();
}
//
//////////////////////////////////////////////////////////////////////
//
namespace
{
    // Error return values for the functions below
    const int ERROR_NO_EXPON = 999;
    const size_t ERROR_NO_INDEX = 999999;
    const size_t ERROR_AM_INDEX = 999998;

    // This function returns the wavelength exponent needed to convert a per wavelength / per frequency
    // quantity to internal (per wavelength) flavor, given the input units, or the error value if the
    // given units are not supported by any of the relevant quantities.
    int waveExponentForSpecificQuantity(Units* unitSystem, string unitString)
    {
        // a list of known per wavelength / per frequency quantities and the corresponding exponents
        static const vector<string> specificQuantities({
                           "wavelengthmonluminosity", "wavelengthfluxdensity", "wavelengthsurfacebrightness",
                           "neutralmonluminosity", "neutralfluxdensity", "neutralsurfacebrightness",
                           "frequencymonluminosity", "frequencyfluxdensity", "frequencysurfacebrightness"});
        static const vector<int> specificExponents({0,0,0, -1,-1,-1, -2,-2,-2});

        // loop over the list
        for (size_t q=0; q!=specificQuantities.size(); ++q)
        {
            // if this quantity supports the given unit, return the corresponding exponent
            if (unitSystem->has(specificQuantities[q], unitString)) return specificExponents[q];
        }
        return ERROR_NO_EXPON;
    }
}

//////////////////////////////////////////////////////////////////////
size_t HDF5InFile::indexForName(std::string name) const
{
    size_t result = ERROR_NO_INDEX;
    size_t index = 0;
    for (const ColumnInfo& col : _colv)
    {
        if (col.title == name)
        {
            if (result != ERROR_NO_INDEX) return ERROR_AM_INDEX;
            result = index;
        }
        index++;
    }
    return result;
}

//////////////////////////////////////////////////////////////////////

size_t HDF5InFile::waveIndexForSpecificQuantity() const
{
    size_t index = 0;
    for (const ColumnInfo& col :_colv)
    {
        if (col.description == "wavelength") return index;
        index++;
    }
    return ERROR_NO_INDEX;
}

//////////////////////////////////////////////////////////////////////
//
void HDF5InFile::useColumns(string columns)
{
    // empty columns string behaves as if we were never called at all
    columns = StringUtils::squeeze(columns);
    if (columns.empty()) return;

    // verify that program columns have not yet been added
    if (_hasProgInfo)
        throw FATALERROR("Program columns were declared before requesting column remapping");

    // verify that file contains column info
    if (!_hasFileInfo)
        throw FATALERROR("Requesting logical columns but there is no column info in file header");

    // establish the logical column info list
    vector<ColumnInfo> newcolv;
    for (string name : StringUtils::split(columns, ","))
    {
        string sname = StringUtils::squeeze(name);
        size_t index = indexForName(sname);
        if (index == ERROR_NO_INDEX)
            throw FATALERROR("No column description in file header matches logical name '" + sname + "'");
        if (index == ERROR_AM_INDEX)
            throw FATALERROR("Multiple column descriptions in file header match logical name '" + sname + "'");

        newcolv.emplace_back(_colv[index]);
    }

    // replace the column info list
    _colv = newcolv;
}

//////////////////////////////////////////////////////////////////////

void HDF5InFile::addColumn(string description, string quantity, string defaultUnit)
{
    _hasProgInfo = true;
    
    size_t index = indexForName(description);

    // get a writable reference to the column record being handled, and increment the program column index
    ColumnInfo& col = _colv[index];
    _numLogCols++;

    // store the programmatically provided information in the record (unit is already stored)
    col.description = description;
    col.quantity = quantity;

    // verify units and determine conversion factor for this column
    if (col.quantity.empty())       // dimensionless quantity
    {
        if (!col.unit.empty() && col.unit != "1")
            throw FATALERROR("Invalid units for dimensionless quantity in column " + std::to_string(_numLogCols));
        col.unit = "1";
    }
    else if (col.quantity == "specific")    // arbitrarily scaled value per wavelength or per frequency
    {
        col.waveExponent = waveExponentForSpecificQuantity(_units, col.unit);
        if (col.waveExponent == ERROR_NO_EXPON)
            throw FATALERROR("Invalid units for specific quantity in column " + std::to_string(_numLogCols));
        if (col.waveExponent)
        {
            col.waveIndex = waveIndexForSpecificQuantity();
            if (col.waveIndex == ERROR_NO_INDEX)
                throw FATALERROR("No preceding wavelength column for specific quantity in column "
                                                       + std::to_string(_numLogCols));
        }
    }
    else
    {
        if (!_units->has(col.quantity, col.unit))
        {
            throw FATALERROR("Invalid units for quantity in column " + std::to_string(_numLogCols));
        }
        col.convFactor = _units->in(col.quantity, col.unit, 1.);
    }
    // add the physical to logical column mapping for this column
    if (_logColIndices.size() < col.physColIndex) _logColIndices.resize(col.physColIndex, ERROR_NO_INDEX);
    if (_logColIndices[col.physColIndex-1] != ERROR_NO_INDEX)
        throw FATALERROR("Multiple logical columns (" + std::to_string(_logColIndices[col.physColIndex-1]+1)
                            + "," + std::to_string(_numLogCols) + ") map to the same physical column ("
                         + std::to_string(col.physColIndex) + ")");
    _logColIndices[col.physColIndex-1] = _numLogCols-1;

    // log column information
    string message = "  Column " + std::to_string(_numLogCols) + ": " + col.description + " (" + col.unit + ")";
    if (!col.title.empty())
    {
        message += " <-- ";
        if (col.physColIndex != _numLogCols) message += "column " + std::to_string(col.physColIndex) + ": ";
        message += col.title;
    }
    _log->info(message);
}

//////////////////////////////////////////////////////////////////////
//
bool HDF5InFile::readRow(Array& values)
{
   readData();
   if (!_hasProgInfo) throw FATALERROR("No columns were declared for column text file");
    
   if(_numRows == _currentRowIndex)
       return false;
   // read new line until it is non-empty and non-comment
       
   // resize result array if needed (we don't need it to be cleared)
   if (values.size() != _numLogCols) values.resize(_numLogCols);
    for (size_t i : _logColIndices)         // i: zero-based logical index
        {
            // read the value as floating point
            double value = _data[i][_currentRowIndex];
            
            // if mapped to a logical column, convert from input units to internal units, and store the result
            if (i != ERROR_NO_INDEX)
            {
                const ColumnInfo& col = _colv[i];
                values[i] = value * (col.waveExponent ? pow(values[col.waveIndex],col.waveExponent)
                                                        : col.convFactor);
            }
        }
    _currentRowIndex++;
    return true;
}
//
//////////////////////////////////////////////////////////////////////
//
bool HDF5InFile::readNonLeaf(int& nx, int& ny, int& nz)
{
    throw FATALERROR("readNonLeaf not implemented");
//    string line;
//
//    while (true)
//    {
//        int c = _in.peek();
//
//        // skip comments line
//        if (c=='#')
//        {
//            getline(_in,line);
//        }
//
//        // process nonleaf line
//        else if (c=='!')
//        {
//            _in.get();              // skip exclamation mark
//            getline(_in,line);
//
//            // convert nx,ny,nz values from line and store them in output arguments
//            std::stringstream linestream(line);
//            linestream >> nx >> ny >> nz;
//            if (linestream.fail())
//                throw FATALERROR("Nonleaf subdivision specifiers are missing or not formatted as integers");
//
//            return true;
//        }
//
//        // eat leading white space and empty lines
//        else if (c==' ' || c=='\t' || c=='\n' || c=='\r')
//        {
//            _in.get();
//        }
//
//        // signal not a nonleaf line
//        else
//        {
//            return false;
//        }
//    }
}
//
//////////////////////////////////////////////////////////////////////
//
vector<Array> HDF5InFile::readAllRows()
{
   vector<Array> rows;
   while (true)
   {
       rows.emplace_back();        // add a default-constructed array to the vector
       if (!readRow(rows.back()))  // read next line's values into that array
       {
           rows.pop_back();        // at the end, remove the extraneous array
           break;
       }
   }
   return rows;
}
//
//////////////////////////////////////////////////////////////////////
//
vector<Array> HDF5InFile::readAllColumns()
{
    
   // read the remainder of the file into rows
   const vector<Array>& rows = readAllRows();
   size_t nrows = rows.size();
   size_t ncols = _colv.size();

   // transpose the result into columns
   vector<Array> columns(ncols, Array(nrows));
   for (size_t c=0; c!=ncols; ++c)
       for (size_t r=0; r!=nrows; ++r)
           columns[c][r] = rows[r][c];

   return columns;
}
//
//////////////////////////////////////////////////////////////////////


void HDF5InFile::readData()
{
    if(!_data.empty())
        return;
    else
    {
        size_t expectedNumRows = 0;
        _data.resize(_colv.size());
        for (const ColumnInfo& col :_colv)
        {
            size_t myIndex = _logColIndices[col.physColIndex-1];
            _data[myIndex] = H5Easy::load<std::vector<double> >(*_in,col.title);
            size_t myNumRows = _data[myIndex].size();
            if(expectedNumRows==0)
                expectedNumRows = myNumRows;
            else
            {
                if(expectedNumRows != myNumRows)
                {
                    throw FATALERROR("The number of rows in each HDF5 dataset needs to be the same!");
                }
            }
        }
        _numRows = expectedNumRows;
    }
    
}
