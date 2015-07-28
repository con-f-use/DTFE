/*
 *  Copyright (c) 2011       Marius Cautun
 *
 *                           Kapteyn Astronomical Institute
 *                           University of Groningen, the Netherlands
 *
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */



/* This file defines a reader for binary Gadget snapshot files and also for HDF5 files (need to have the HDF5 library for that). The reader given here can read all the particles and properties stored in the Gadget snapshot file by giving the appropriate options to the program. See the documentation for the options that control this behavior. */



#define SWAP_HEADER_ENDIANNESS(x1,x2,x3,x4) { if( x1 ) {BYTESWAP( x2 ); BYTESWAP( x3 ); x4.swapBytes();} }
#define SWAP_ENDIANNESS(x1,x2,x3)           { if( x1 ) {BYTESWAP( x2 ); BYTESWAP( x3 );} }
// do ... while(0) in the following is hack to emulate regular expression behavior in conditional clauses
#define READ_DELIMETER do \
{\
    inputFile.seekg( offset, std::ios::cur );\
    inputFile.read( reinterpret_cast<char *>(&buffer1), sizeof(buffer1) );\
    if( swapEndian ) BYTESWAP( buffer1 ); \
} while(0);
#define DELIMETER_CONSISTANCY_CHECK(field) do \
{\
    inputFile.read( reinterpret_cast<char *>(&buffer2), sizeof(buffer2) ); \
    if( swapEndian ) BYTESWAP( buffer2 ); \
    if ( buffer1!=buffer2 ) { \
        message << "\nbefore: " << buffer1 << " after: " << buffer2 << "\n" << MESSAGE::Flush;\
        throwError( "The integers before and after the particle " field " data block in the GADGET file '" + fileName + "' did not match. The GADGET snapshot file is corrupt." ); \
    } \
} while(0);

// Header structure for reading Gadget snapshots
struct Gadget_header
{
    int      npart[6];
    double   mass[6];
    double   time;
    double   redshift;
    int      flag_sfr;
    int      flag_feedback;
    int      npartTotal[6];
    int      flag_cooling;
    int      num_files;
    double   BoxSize;
    double   Omega0;
    double   OmegaLambda;
    double   HubbleParam;
    char     fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8];  /* fills to 256 Bytes */


    // return the file name for a Gadget snapshot saved in single or multiple files - note that the name must contain a '%i' or '%s' character
    std::string filename(std::string fileRoot, int const fileNumber, bool checkFileExists=true )
    {
        char buf[500];
        sprintf( buf, fileRoot.c_str(), fileNumber );
        std::string fileName( buf );
        if ( not bfs::exists(fileName) and checkFileExists )
            throwError( "The program could not open the input GADGET snapshot file/files: '" + fileName + "'. It cannot find the file/files." );
        return fileName;
    }

    // Function that prints the Gadget header.
    void print()
    {
        std::cout << "\nThe header of the Gadget file contains the following info:\n"
            << "npart[6]     =  " << npart[0] << "  " << npart[1] << "  " << npart[2] << "  " << npart[3] << "  " <<  npart[4] << "  " <<  npart[5] << "\n"
            << "mass[6]      =  " << mass[0] << "  " << mass[1] << "  " << mass[2] << "  " << mass[3] << "  " << mass[4] << "  " << mass[5] << "\n"
            << "time         =  " << time << "\n"
            << "redshift     =  " << redshift << "\n"
            << "flag_sfr     =  " << flag_sfr << "\n"
            << "flag_feedback=  " << flag_feedback << "\n"
            << "npartTotal[6]=  " << npartTotal[0] << "  " << npartTotal[1] << "  " << npartTotal[2] << "  " << npartTotal[3] << "  " << npartTotal[4] << "  " << npartTotal[5] << "  " << "\n"
            << "flag_cooling =  " << flag_cooling << "\n"
            << "num_files    =  " << num_files << "\n"
            << "BoxSize      =  " << BoxSize << "\n"
            << "Omega0       =  " << Omega0 << "\n"
            << "OmegaLambda  =  " << OmegaLambda << "\n"
            << "h            =  " << HubbleParam << "\n\n";
    }

    // Swap endianness
    void swapBytes()
    {
        ByteSwapArray( npart, 6 );
        ByteSwapArray( mass, 6 );
        BYTESWAP( time );
        BYTESWAP( redshift );
        BYTESWAP( flag_sfr );
        BYTESWAP( flag_feedback );
        ByteSwapArray( npartTotal, 6 );
        BYTESWAP( flag_cooling );
        BYTESWAP( num_files );
        BYTESWAP( BoxSize );
        BYTESWAP( Omega0 );
        BYTESWAP( OmegaLambda );
        BYTESWAP( HubbleParam );
    }

    // Checks for the type of the Gadget file -> can detected Gadget file type 1 & 2. Returns true if it could identify the gadget file type.
    bool detectSnapshotType(int const bufferValue,
                            int *gadgetFileType,
                            bool *swapEndian)
    {
        int buffer1 = bufferValue;
        *swapEndian = false;

        if ( buffer1 == 8 )             // gadget file format 2
            *gadgetFileType = 2;
        else if ( buffer1 == 256 )      // gadget file format 1
            *gadgetFileType = 1;
        else                            // check for swapped endianness
        {
            BYTESWAP( buffer1 );
            *swapEndian = true;
            if ( buffer1 == 8 )         // gadget file format 2
                *gadgetFileType = 2;
            else if ( buffer1 == 256 )  // gadget file format 1
                *gadgetFileType = 1;
            else                        // could not detect the file type
                return false;
        }
        return true;
    }
};




// function that counts the number of particles for snapshots saved as multiple files
void countGadgetParticleNumber(std::string filenameRoot,
                               int const noFiles,
                               int const gadgetFileType,
                               bool const swapEndian,
                               int const verboseLevel,
                               size_t numberTotalParticles[]);

// function that reads the data in a single gadget file - this function is called by the 'readGadgetFile_single' and 'readGadgetFile_multiple' to do the actual data reading
void readGadgetData(std::string fileName,
                    Read_data<float> *readData,
                    User_options &userOptions,
                    int const gadgetFileType,
                    bool const swapEndian,
                    size_t *numberParticlesRead);




/*! This functions reads the gadget header for one of the input files and uses that to set all the properties need for reading the input data:
        - finding the number of files
        - finding the header type and endianess of the data
        - setting the box length in the header as dimensions of the box
        - finding the number of particles to be read from the input files
        - reserving memory for reading the data
 */
void initializeGadget(std::string filename,
                      Read_data<float> *readData,
                      User_options *userOptions,
                      Gadget_header *gadgetHeader,
                      int *gadgetFileType,
                      bool *swapEndian,
                      size_t *noParticles)
{
    MESSAGE::Message message( userOptions->verboseLevel );
    std::string fileName = filename;
    bool singleFile = true;
    if ( not bfs::exists(fileName) ) //if this is true, than the input is in several files
    {
        fileName = gadgetHeader->filename( filename, 0 );
        singleFile = false;
    }


    // open the first binary file for reading and read some of the overall characteristics
    std::fstream inputFile;
    openInputBinaryFile( inputFile, fileName );


    // detect the Gadget file type -> gadget file type 1 or 2
    int buffer1, buffer2;       // variables to read the buffer before and after each gadget data block
    inputFile.read( reinterpret_cast<char *>(&buffer1), sizeof(buffer1) );
    bool validFile = gadgetHeader->detectSnapshotType( buffer1, gadgetFileType, swapEndian );
    if ( not validFile )
        throwError( "Unknown file type for the input Gadget snapshot. Tried Gadget snapshots type 1 and 2 as well as changing endianness, but none worked. Check that you inserted the correct input file." );
    if ( *swapEndian )
        message << "Detected that the input data file has a different endianness than the current system. The program will automatically change endianness for the data!" << MESSAGE::Flush;
    int offset = (*gadgetFileType)==2 ? 16 : 0;      // keep track if file type 2 to have a 16 bytes offset every time reading a new data block


    // now read the actual values of the gadget header
    inputFile.seekg( offset, std::ios::beg );
    inputFile.read( reinterpret_cast<char *>(&buffer1), sizeof(buffer1) );
    inputFile.read( reinterpret_cast<char *>(gadgetHeader), sizeof(*gadgetHeader) );
    inputFile.read( reinterpret_cast<char *>(&buffer2), sizeof(buffer2) );
    inputFile.close();
    SWAP_HEADER_ENDIANNESS( *swapEndian, buffer1, buffer2, (*gadgetHeader) ); //swap endianness if that is the case
    if ( buffer1!=buffer2 or buffer1!=256 )
        throwError( "The was an error while reading the header of the GADGET snapshot file. The integers before and after the header do not match the value 256. The GADGET snapshot file is corrupt." );


    // check if to set the box coordinates from the information in the header
    if ( not userOptions->userGivenBoxCoordinates )
    {
        for (size_t i=0; i<NO_DIM; ++i)
        {
            userOptions->boxCoordinates[2*i] = 0.;                    // this is the left extension of the full box
            userOptions->boxCoordinates[2*i+1] = gadgetHeader->BoxSize;// right extension of the full box
        }
    }
    else
        message << "The box coordinates were set by the user using the program options. The program will keep this values and will NOT use the box length information from the Gadget file!" << MESSAGE::Flush;
#ifdef WOJTEK
    if ( userOptions->additionalOptions.size()!=0 ) //if inserted an option
        gadgetHeader->num_files = atoi( userOptions->additionalOptions[0].c_str() );
    gadgetHeader->print();
#endif


    // get the total number of particles to be read from the file/files
    size_t numberTotalParticles[6];
    if ( singleFile )    // if single file
    {
        gadgetHeader->num_files = 1;
        for (int i=0; i<6; ++i) {
            numberTotalParticles[i] = gadgetHeader->npart[i];
        }
    }
    else                  // for multiple files read the number of particles in each file and keep track of that
        countGadgetParticleNumber( filename, gadgetHeader->num_files, *gadgetFileType, *swapEndian, userOptions->verboseLevel, numberTotalParticles );

    // from that keep only the particles specified by the user (not all the particle species may be of interest)
    *noParticles = 0;   //total number of particles to be read from file
    for (int i=0; i<6; ++i)
    {
        if ( not userOptions->readParticleSpecies[i] )
            numberTotalParticles[i] = 0;
        *noParticles += numberTotalParticles[i];
    }
    message << "The data consists of " << gadgetHeader->npart[0];
    size_t noAllPart = gadgetHeader->npart[0];
    for (int i=1; i<6; ++i) {
        noAllPart += gadgetHeader->npart[i];
        message << " + " << gadgetHeader->npart[i];
    }
    message << " = " << noAllPart << " particles listed by type. Reading a total of " << *noParticles << " of these particles as follows: " << numberTotalParticles[0];
    for (int i=1; i<6; ++i)
        message << " + " << numberTotalParticles[i];
    message << ".\n" << MESSAGE::Flush;



    // allocate memory for the particle data
    message << "Allocating memory for " << MESSAGE::Flush;
    if ( userOptions->readParticleData[0] ) {
        message << "positions (" << *noParticles << ")... " << MESSAGE::Flush;
        readData->position( *noParticles );  // particle positions
    }
    if ( userOptions->readParticleData[1] ) {
        message << "weights (" << *noParticles << ")... " << MESSAGE::Flush;
        readData->weight( *noParticles );    // particle weights (weights = particle masses from the snapshot file)
    }
#ifdef VELOCITY
    if ( userOptions->readParticleData[2] ) {
        message << "velocity (" << *noParticles << ")... " << MESSAGE::Flush;
        readData->velocity( *noParticles );  // particle velocities
    }
#endif
#ifdef SCALAR
    int noScalars = 0; // Number of scalar fields
    for (size_t i=3; i<userOptions->readParticleData.size(); ++i)
        if ( userOptions->readParticleData[i] )
            noScalars += 1;
    if ( noScalars>0 ) {
        message << "scalars... " << MESSAGE::Flush;
        readData->scalar( *noParticles );  // particle scalar quantity
    }
    message << "Done.\n";
#endif

    // check that the 'readParticleData' and 'readParticleSpecies' options make sense
    if ( not userOptions->readParticleData[0] )
        throwError( "The program needs the particle position information to be able to interpolate the fields on a grid. Please add '1' to the integer number giving the data blocks to be read from the input Gadget snapshot." );
    if ( not userOptions->readParticleData[1] )
    {
        MESSAGE::Warning warning( userOptions->verboseLevel );
        warning << "You selected not to read the Gadget particle masses. This means that the program will treat all particles as having the same weight (mass)." << MESSAGE::EndWarning;
    }
    if ( (*noParticles)<=0 )
        throwError( "Please select again the particle species that you would like to read from the Gadget file. There are no particles in the current selection!" );
}



/*! This function reads a Gadget snapshots saved in a single or multiple files.

NOTE: This function is a not the best choice to understand how to read input data for users unfamiliar with the Gadget snapshot files. */
void readGadgetFile(std::string filename,
                    Read_data<float> *readData,
                    User_options *userOptions)
{
    int gadgetFileType;     // stores the type of the gadget file (1 or 2)
    bool swapEndian;        // true if need to swap the endianess of the data
    size_t noParticles;     // total number of particles to be read from the file
    Gadget_header gadgetHeader; // the gadget header for the first file


    // get the input gadget data characteristics: file type, endianess, particle number and number of files. The following function also reserves memory for the positions, weights and velocity data (if requested so).
    initializeGadget( filename, readData, userOptions, &gadgetHeader, &gadgetFileType, &swapEndian, &noParticles );


    MESSAGE::Message message( userOptions->verboseLevel );
    std::string fileName = filename;


    // read the data in each file
    size_t numberParticlesRead = 0;   // the total number of particles read after each file

    // iterate over all the files and read the data
    for (int i=0; i<gadgetHeader.num_files; ++i )
    {
        fileName = gadgetHeader.filename( filename, i );
        message << "Reading GADGET snapshot file '" << fileName << "' which is file " << i+1 << " of " << gadgetHeader.num_files << " files...\n" << MESSAGE::Flush;

        // call the function that reads the data
        readGadgetData( fileName, readData, *userOptions, gadgetFileType, swapEndian, &numberParticlesRead );
    }


    // check if the data needs to be swapped
    if ( swapEndian )
    {
        if ( userOptions->readParticleData[0] )
            ByteSwapArray( readData->position(), NO_DIM*noParticles );
        if ( userOptions->readParticleData[1] )
            ByteSwapArray( readData->weight(), noParticles );
#ifdef VELOCITY
        if ( userOptions->readParticleData[2] )
            ByteSwapArray( readData->velocity(), NO_DIM*noParticles );
#endif
    }

    return;
}




/*! Function that counts how many particles are in a Gadget file snapshot that was saved in multiple files.
*/
void countGadgetParticleNumber(std::string filenameRoot,
                               int const noFiles,
                               int const gadgetFileType,
                               bool const swapEndian,
                               int const verboseLevel,
                               size_t numberTotalParticles[])
{
    for (int i=0; i<6; ++i)
        numberTotalParticles[i] = 0;
    int offset = gadgetFileType==2 ? 16 : 0;

    // loop over all the files and count the number of particles
    for (int i=0; i<noFiles; ++i)
    {
        Gadget_header gadgetHeader;
        std::string fileName = gadgetHeader.filename( filenameRoot, i ); //get the filename for file i

        std::fstream inputFile;
        openInputBinaryFile( inputFile, fileName );

        // read the header
        int buffer1, buffer2;
        inputFile.seekg( offset, std::ios::beg );
        inputFile.read( reinterpret_cast<char *>(&buffer1), sizeof(buffer1) );
        inputFile.read( reinterpret_cast<char *>(&gadgetHeader), sizeof(gadgetHeader) );
        inputFile.read( reinterpret_cast<char *>(&buffer2), sizeof(buffer2) );
        inputFile.close();
        SWAP_HEADER_ENDIANNESS( swapEndian, buffer1, buffer2, gadgetHeader ); //swap endianness if that is the case
        if ( buffer1!=buffer2 or buffer1!=256 )
            throwError( "The was an error while reading the header of the GADGET snapshot file '" + fileName + "'. The integers before and after the header do not match the value 256. The GADGET snapshot file is corrupt." );

        // add the particle numbers
        for (int j=0; j<6; ++j)
            numberTotalParticles[j] += gadgetHeader.npart[j];
    }

    MESSAGE::Message message( verboseLevel );
    message << "The data is in " << noFiles << " files and contains the following number of particles: " << numberTotalParticles[0] << " + "  << numberTotalParticles[1] << " + "  << numberTotalParticles[2] << " + "  << numberTotalParticles[3] << " + "  << numberTotalParticles[4] << " + "  << numberTotalParticles[5] << " .\n" << MESSAGE::Flush;
}



/*! This function reads the gadget data from a single file. This function should be called from within a loop to read a gadget snapshot saved in multiple files.
*/
void readGadgetData(std::string fileName,
                    Read_data<float> *readData,
                    User_options &userOptions,
                    int const gadgetFileType,
                    bool const swapEndian,
                    size_t *numberParticlesRead)
{
    MESSAGE::Message message( userOptions.verboseLevel );
    int offset = gadgetFileType == 2 ? 16 : 0;


    // open the file
    std::fstream inputFile;
    openInputBinaryFile( inputFile, fileName );


    // read block #1 header (block #1 from gadget's manual)
    int buffer1, buffer2;
    Gadget_header tempHeader;
    READ_DELIMETER;
    inputFile.read( reinterpret_cast<char *>(&tempHeader), sizeof(tempHeader) );
    DELIMETER_CONSISTANCY_CHECK("header");
    if ( buffer1!=256 )
        throwError( "The integers before and after the header do not match the value 256. The GADGET snapshot file is corrupt." );


    // read block #2 positions
    READ_DELIMETER;
    if ( userOptions.readParticleData[0] )
    {
        float *positions = readData->position();               // returns a pointer to the particle positions array
        size_t dataOffset = (*numberParticlesRead) * NO_DIM;    // the offset in the positions array from where to start reading the new positions
        message << "\treading positions of the particles starting at " << dataOffset << "... " << MESSAGE::Flush;

        // loop over each species and reads its data if the user requested so
        for (int i=0; i<6; ++i)
        {
            // check if to skip this data block
            if ( tempHeader.npart[i]==0 )       // skip since no particle of this species is present
            {
                ;
            }
            else if ( not userOptions.readParticleSpecies[i] )    // skip this data block since don't need this data
            {
                inputFile.seekg( tempHeader.npart[i] * sizeof(float) * NO_DIM, std::ios::cur );
            }
            // read this data block
            else
            {
                inputFile.read( reinterpret_cast<char *>( &(positions[dataOffset]) ), tempHeader.npart[i] * sizeof(float) * NO_DIM );
                dataOffset += tempHeader.npart[i] * NO_DIM;
            }
        }

        // check for consistency in the input file
        message << "Done.";
    }
    else    //skip the positions if not interested in them
    {
        message << "\n\t(skipping positions)" << MESSAGE::Flush;
        inputFile.seekg( buffer1, std::ios::cur );
    }
    DELIMETER_CONSISTANCY_CHECK("position");


    // read block #3 velocities
    READ_DELIMETER;
#ifdef VELOCITY
    if ( userOptions.readParticleData[2] )
    {
        float *velocities = readData->velocity();              // returns a pointer to the particle velocity array
        size_t dataOffset = (*numberParticlesRead) * NO_DIM;    // the offset in the velocities array from where to start reading the new velocities
        message << "\n\treading velocities of the particles... " << MESSAGE::Flush;

        // loop over each species and reads its data if the user requested so
        for (int i=0; i<6; ++i)
        {
            // check if to skip this data block
            if ( tempHeader.npart[i]==0 )       // skip since no particle of this species is present
            {
                ;
            }
            // skip this data block since don't need this data
            else if ( not userOptions.readParticleSpecies[i] )
            {
                inputFile.seekg( tempHeader.npart[i] * sizeof(float) * NO_DIM, std::ios::cur );
            }
            // read this data block
            else
            {
                inputFile.read( reinterpret_cast<char *>( &(velocities[dataOffset]) ), tempHeader.npart[i] * sizeof(float) * NO_DIM );
                dataOffset += tempHeader.npart[i] * NO_DIM;
            }
        }
        message << "Done.";
    }
    else
    {   //skip the velocities if not interested in them
        message << "\n\tskipping velocities... " << MESSAGE::Flush;
        inputFile.seekg( buffer1, std::ios::cur );
        message << "Done." << MESSAGE::Flush;
    }
#else
    message << "\n\tskipping velocities... " << MESSAGE::Flush;
    inputFile.seekg( buffer1, std::ios::cur );
    message << "Done." << MESSAGE::Flush;
#endif
    DELIMETER_CONSISTANCY_CHECK("velocity");


    // skip block #4 particle IDs
    message << "\n\tskipping ids... " << MESSAGE::Flush;
    READ_DELIMETER;
    inputFile.seekg( buffer1, std::ios::cur );
    message << "Done." << MESSAGE::Flush;
    DELIMETER_CONSISTANCY_CHECK("id");


    // read block #5 masses (or weights)
    bool massBlockPresent = false;
    for (int i=0; i<6; ++i)
        if ( tempHeader.mass[i]==0. and tempHeader.npart[i]!=0 )
            massBlockPresent = true;

    if ( massBlockPresent ) READ_DELIMETER;
    if ( userOptions.readParticleData[1] )
    {
        float *weights = readData->weight();          // returns a pointer to the particle weights array
        size_t dataOffset = (*numberParticlesRead);    // the offset in the masses array from where to start reading the new masses
        size_t readinthisiterations = 0;
        message << "\n\treading masses of the particles starting at " << *numberParticlesRead << "... \n" << MESSAGE::Flush;

        // loop over each species and reads its data if the user requested so
        for (int i=0; i<6; ++i)
        {
            message << "\t\tType " << i+1 << " (" << tempHeader.npart[i] << " particles): " << MESSAGE::Flush;
            // check if to skip this data block
            if ( tempHeader.npart[i]==0 )       // skip since no particle of this species is present
            {
                message << "No Particles.\n" << MESSAGE::Flush;
            }
            else if ( tempHeader.mass[i] <= 0.0 )  // particles have individual masses
            {
                // read masses
                if ( userOptions.readParticleSpecies[i] ) {
                    message << "Have individual masses..." << MESSAGE::Flush;
                    inputFile.read( reinterpret_cast<char *>( &(weights[dataOffset]) ), tempHeader.npart[i] * sizeof(float) );
                    dataOffset += tempHeader.npart[i];
                    message << " complete.\n" << MESSAGE::Flush;
                    readinthisiterations += tempHeader.npart[i];
                }
                else
                {
                    message << "Skipped by user (individual mass).\n" << MESSAGE::Flush;
                    inputFile.seekg( tempHeader.npart[i] * sizeof(float), std::ios::cur );
                    readinthisiterations += tempHeader.npart[i];
                }
            }
            else if ( userOptions.readParticleSpecies[i] )
            {
                message << "All have the same masses of " << tempHeader.mass[i] << MESSAGE::Flush;
                float mass = tempHeader.mass[i];
                if ( swapEndian ) BYTESWAP(mass);       // keep the same endianess for all the data
                for (size_t j=dataOffset; j<dataOffset+tempHeader.npart[i]; ++j) {
                    weights[j] = mass;
                }
                dataOffset += tempHeader.npart[i];
                message << " complete.\n" << MESSAGE::Flush;
            }
            else
            {
                message << "Skipped by user (per particle mass from header: " << tempHeader.mass[i] << ")\n" << MESSAGE::Flush;
            }
        }
        message << "\t...read or skipped " << readinthisiterations << " particles in this go. Done." << MESSAGE::Flush;
    }
    if ( massBlockPresent ) DELIMETER_CONSISTANCY_CHECK("mass");


#ifdef SCALAR
    // prepare to read scalar properties
    int noScalars = 0; // Number of scalar fields
    for (size_t i=3; i<userOptions.readParticleData.size(); ++i)
        if ( userOptions.readParticleData[i] )
            noScalars += 1;
    if ( noScalars>0 ) {

        size_t noScalarsRead = 0;
        float *scalar = readData->scalar();          // returns a pointer to the particle scalar properties array


        // read block #6 the internal energies U
        if( userOptions.readParticleData[3] and ( inputFile.eof() or tempHeader.npart[0]<1 ) )
            message << "\t no particles with internal energy found in '" << fileName << "'";

        if( !inputFile.eof() ) // from here on the file might be at an end
        {
            READ_DELIMETER;
            if ( buffer1/sizeof(float) == tempHeader.npart[0] ) { // it can only be an int. energy block if the number of floats in the block is equal to the number of gas particles
                if ( userOptions.readParticleData[3] )
                {
                    float *weights = readData->weight();         // returns a pointer to the particle weights array
                    message << "\t reading internal energy of the particles... " << MESSAGE::Flush;

                    // read the energies
                    float *tempData = new float[ tempHeader.npart[0] ];
                    inputFile.read( reinterpret_cast<char *>( tempData ), tempHeader.npart[0] * sizeof(float) );

                    float mean = 0;
                    for (size_t i; i<size_t(tempHeader.npart[0]); ++i)
                    {
                        size_t index1 = (*numberParticlesRead) + i;
                        size_t index2 = index1 * NO_SCALARS + noScalarsRead;
                        scalar[index2] = weights[index1] * tempData[i]; // U is given per unit mass in Gadget
                        mean += scalar[index2];
                    }
                    mean /= tempHeader.npart[0];
                    // or read file on the fly (more memory efficient but slower): //for (size_t i=(*numberParticlesRead); i<size_t(tempHeader.npart[0])+(*numberParticlesRead); ++i ) { float tempData; inputFile.read( reinterpret_cast<char *>( &tempData ), sizeof(tempData) ); scalar[i*NO_SCALARS+noScalarsRead] = weights[i]*tempData; }

                    delete[] tempData;
                    noScalarsRead++;
                    message << "\r\t reading internal energy of the particles (mean energy: " << mean << ")... Done." << MESSAGE::Flush;
                    DELIMETER_CONSISTANCY_CHECK("internal energy U");
                }
                else
                    // no skipping here, because: only the internal energy is given per unit mass and must be treated speciallyif the user does not read U we can treat (and skip) it below. Normally U corresponds to 8 in the "--input" options, the block after U corresponds to 16. If U is present in the snapshot and the user has not specified 8 in the "--input" the U block will be treated as an unknown gas block so as 16 and the block after U will correspond to 32. This way if there is no U block but a gas block, we avoid for the code to falsely recognize it as U block and multiply by the mass.
                    inputFile.seekg( -sizeof(buffer1), std::ios::cur ); // revert reading the buffer so code afterwards can continue normally
            }
            else
                inputFile.seekg( -sizeof(buffer1), std::ios::cur ); // see above
        } else if ( userOptions.readParticleData[3] )
            throwError("Expected internal energy block U was not present in '", fileName,"'!");

        // read other blocks
        // A total of 12 blocks of data is mentioned in gadget's manual, which leaves additional 6 blocks
        // (Note: Up to (currently) 40 blocks may be present in the full physics version of gadget which are not yet read out by the code)
        for(int i=4; i<userOptions.readParticleData.size(); ++i)
      //for(int i=4; 1; ++i) if all blocks shell be process but there are only 10 entries in readParticleData
        {
            if( !inputFile.eof() ) {
                READ_DELIMETER;
                size_t readBytes = 0;
                for (int j=0; j<6; ++j) readBytes += tempHeader.npart[j] * sizeof(float);
                if( userOptions.readParticleData[i] )
                {
                    message << "\n\t reading additional block " << i;
                    // 1D property of all particles (e.g. potential, timestep)
                    if ( buffer1 == readBytes )
                        message << " of the particles... " << MESSAGE::Flush;
                    // 1D property of gas particles (e.g. metalicity, density, sfr, hsml, dAdt)
                    else if ( buffer1 == tempHeader.npart[0] * sizeof(float) )
                        message << " of gas particles... " << MESSAGE::Flush;
                    // 3D property of all particles
                    else if ( buffer1 == readBytes * NO_DIM )
                    {
                        message << " of the particles (3D)... " << MESSAGE::Flush;
                        throwError("The program does not know how to read 3D scalars yet!");
                    }
                    // none of the above
                    else
                        throwError("Number of particles in header does not match number in '", fileName, "' (", buffer1/sizeof(float),").");

                    float *tempData = new float[ buffer1/sizeof(float) ];
                    inputFile.read( reinterpret_cast<char *>( tempData ), buffer1 );
                    for (size_t k; k<size_t(buffer1/sizeof(float)); ++k)
                    {
                        size_t index1 = ((*numberParticlesRead) + k) * NO_SCALARS + noScalarsRead;
                        scalar[index1] = tempData[k];
                    }
                    delete[] tempData;
                    noScalarsRead++;
                    message << "Done.\n" << MESSAGE::Flush;
                }
                else
                {
                    message << "\t (skipping unk. block #" << i-3 <<")" << MESSAGE::Flush;
                    if(i==5) message << "\n";
                    inputFile.seekg( buffer1, std::ios::cur );
                }
                DELIMETER_CONSISTANCY_CHECK("unknown block");
            }
            else if( userOptions.readParticleData[i] )
                throwError("Expected unknown data block #", i-3, " was not present in '", fileName,"'!");
            else
                break;
        }

    }
#endif
    message << "\n";

    inputFile.close();
    for (int i=0; i<6; ++i)     // update the number of read particles
        if ( userOptions.readParticleSpecies[i] )
            (*numberParticlesRead) += tempHeader.npart[i];
}






// The following are the functions used to read an HDF5 gadget file
#ifdef HDF5
#include <H5Cpp.h>
using namespace H5;


/*! Function to check if attribute exists. */
extern "C"
{
    bool doesAttributeExist(hid_t obj_id, const char* name)
    {
        return( H5Aexists( obj_id, name ) > 0 ? true : false );
    }
}




/*! Reads some of the entries for the Gadget header from an HDF5 file.
NOTE: it does not read all the entries, it only reads the particle number in the file, the mass array, the box length and the number of files per snapshot.
*/
void HDF5_readGadgetHeader(std::string filename,
                           Gadget_header *gadgetHeader)
{
    // the name of the HDF5 file
    const H5std_string FILE_NAME( filename );

    // open the HDF5 file and the header group
    H5File *file = new H5File( FILE_NAME, H5F_ACC_RDONLY, H5P_DEFAULT );
    Group *group = new Group( file->openGroup("/Header") );


    // start reading one header attribute at a time
    std::string name( "NumPart_ThisFile" );
    if ( doesAttributeExist( group->getId(), name.c_str() ) )
        group->openAttribute( name.c_str() ).read( PredType::NATIVE_INT, gadgetHeader->npart );
    else throwError( "No '" + name + "' attribute found in the HDF5 file '" + filename + "'. Cannot continue with the program!" );

    name = "MassTable";
    if ( doesAttributeExist( group->getId(), name.c_str() ) )
        group->openAttribute( name.c_str() ).read( PredType::NATIVE_DOUBLE, gadgetHeader->mass );
    else throwError( "No '" + name + "' attribute found in the HDF5 file '" + filename + "'. Cannot continue with the program!" );

    name = "NumFilesPerSnapshot";
    if ( doesAttributeExist( group->getId(), name.c_str() ) )
        group->openAttribute( name.c_str() ).read( PredType::NATIVE_INT, &(gadgetHeader->num_files) );
    else throwError( "No '" + name + "' attribute found in the HDF5 file '" + filename + "'. Cannot continue with the program!" );

    name = "BoxSize";
    if ( doesAttributeExist( group->getId(), name.c_str() ) )
        group->openAttribute( name.c_str() ).read( PredType::NATIVE_DOUBLE, &(gadgetHeader->BoxSize) );

    name = "NumPart_Total";
    if ( doesAttributeExist( group->getId(), name.c_str() ) )
        group->openAttribute( name.c_str() ).read( PredType::NATIVE_UINT, gadgetHeader->npartTotal );

    name = "Redshift";
    if ( doesAttributeExist( group->getId(), name.c_str() ) )
        group->openAttribute( name.c_str() ).read( PredType::NATIVE_DOUBLE, &(gadgetHeader->redshift) );

    name = "Time";
    if ( doesAttributeExist( group->getId(), name.c_str() ) )
        group->openAttribute( name.c_str() ).read( PredType::NATIVE_DOUBLE, &(gadgetHeader->time) );
    name = "Time_GYR";
    if ( doesAttributeExist( group->getId(), name.c_str() ) )
        group->openAttribute( name.c_str() ).read( PredType::NATIVE_DOUBLE, &(gadgetHeader->time) );

    name = "Omega0";
    if ( doesAttributeExist( group->getId(), name.c_str() ) )
        group->openAttribute( name.c_str() ).read( PredType::NATIVE_DOUBLE, &(gadgetHeader->Omega0) );

    name = "OmegaLambda";
    if ( doesAttributeExist( group->getId(), name.c_str() ) )
        group->openAttribute( name.c_str() ).read( PredType::NATIVE_DOUBLE, &(gadgetHeader->OmegaLambda) );

    name = "HubbleParam";
    if ( doesAttributeExist( group->getId(), name.c_str() ) )
        group->openAttribute( name.c_str() ).read( PredType::NATIVE_DOUBLE, &(gadgetHeader->HubbleParam) );


    // close the group and file
    delete group;
    delete file;
}



/*! Reads the Gadget data from a single HDF5 file. This function should be called to do the actual data reading for snapshots saved in single or multiple files.
*/
void HDF5_readGadgetData(std::string filename,
                         Read_data<float> *readData,
                         User_options &userOptions,
                         int const fileIndex,
                         size_t *numberParticlesRead)
{
    MESSAGE::Message message( userOptions.verboseLevel );

    // get the header associated to this file
    Gadget_header gadgetHeader;
    HDF5_readGadgetHeader( filename, &gadgetHeader );

    //select only the species of interest
    for (int i=0; i<6; ++i)
        if ( not userOptions.readParticleSpecies[i] )
            gadgetHeader.npart[i] = 0;


    // open the HDF5 file
    const H5std_string FILE_NAME( filename );
    H5File *file = new H5File( FILE_NAME, H5F_ACC_RDONLY, H5P_DEFAULT );
    Group *group;


    // read the coordinates
    if ( userOptions.readParticleData[0] )
    {
        float *positions = readData->position();               // returns a pointer to the particle positions array
        size_t dataOffset = (*numberParticlesRead) * NO_DIM;   // the offset in the positions array from where to start reading the new positions
        message << "\t reading the particles positions ... " << MESSAGE::Flush;
        for(int type=0; type<6; type++)                        // loop over the particle species
        {
            if ( gadgetHeader.npart[type]<=0 ) continue;
            char buf[500];
            sprintf( buf, "/PartType%d", type );
            group = new Group( file->openGroup(buf) );

            // open the data set
            DataSet dataset = group->openDataSet("Coordinates");

            dataset.read( &(positions[dataOffset]), PredType::NATIVE_FLOAT );
            delete group;

            dataOffset += gadgetHeader.npart[type] * NO_DIM;
        }
        message << "Done\n";
    }


    // read the masses (or weights) if different
    if ( userOptions.readParticleData[1] )
    {
        float *weights = readData->weight();          // returns a pointer to the particle weights array
        size_t dataOffset = (*numberParticlesRead);    // the offset in the masses array from where to start reading the new masses
        message << "\t reading the particles masses ... " << MESSAGE::Flush;
        for(int type=0; type<6; type++)                        // loop over the particle species
        {
            if ( gadgetHeader.npart[type]<=0 ) continue;
            if ( gadgetHeader.mass[type]!=0. )          // particle masses given in header
            {
                float mass = gadgetHeader.mass[type];
                for (size_t j=dataOffset; j<dataOffset+gadgetHeader.npart[type]; ++j)
                    weights[j] = mass;
            }
            else                                        // particle masses are variable
            {
                char buf[500];
                sprintf( buf, "/PartType%d", type );
                group = new Group( file->openGroup(buf) );

                DataSet dataset = group->openDataSet("Mass");

                dataset.read( &(weights[dataOffset]), PredType::NATIVE_FLOAT );
                delete group;
            }
            dataOffset += gadgetHeader.npart[type];
        }
        message << "Done\n";
    }


    // read the velocities
    if ( userOptions.readParticleData[2] )
    {
        float *velocities = readData->velocity();              // returns a pointer to the particle velocity array
        size_t dataOffset = (*numberParticlesRead) * NO_DIM;   // the offset in the velocity array from where to start reading the new data
        message << "\t reading the particles velocities ... " << MESSAGE::Flush;
        for(int type=0; type<6; type++)                        // loop over the particle species
        {
            if ( gadgetHeader.npart[type]<=0 ) continue;
            char buf[500];
            sprintf( buf, "/PartType%d", type );
            group = new Group( file->openGroup(buf) );

            // open the data set
            DataSet dataset = group->openDataSet("Velocity");

            dataset.read( &(velocities[dataOffset]), PredType::NATIVE_FLOAT );
            delete group;

            dataOffset += gadgetHeader.npart[type] * NO_DIM;
        }
        message << "Done\n";
    }



    size_t noScalarsRead = 0;
    // read the temperatures
    if ( userOptions.readParticleData[3] and gadgetHeader.npart[0]>0 )
    {
        // the tempartures make sense only for gas particles
        group = new Group( file->openGroup( "/PartType0" ) );
        DataSet dataset = group->openDataSet("Temperature");
        float *tempData = new float[ gadgetHeader.npart[0] ];
        dataset.read( tempData, PredType::NATIVE_FLOAT );
        delete group;

        // get the mass weighted temperature
        float *scalar = readData->scalar();          // returns a pointer to the particle scalar properties array
        float *weights = readData->weight();         // returns a pointer to the particle weights array
        size_t dataOffset = (*numberParticlesRead);  // the offset in the scalar array from where to start reading the new values
        for (size_t i=0; i<size_t(gadgetHeader.npart[0]); ++i)
        {
            size_t index1 = dataOffset + i;
            size_t index2 = index1 * NO_SCALARS + noScalarsRead;
            scalar[index2] = weights[index1] * tempData[i];
        }

        delete[] tempData;
        noScalarsRead += 1;
    }
    // read the HI content
    if ( userOptions.readParticleData[4] and gadgetHeader.npart[0]>0 )
    {
        // the element abundance makes sense only for gas particles
        group = new Group( file->openGroup( "/PartType0/ElementAbundance" ) );
        DataSet dataset = group->openDataSet("Hydrogen");
        float *hydrogenMassFraction = new float[ gadgetHeader.npart[0] ];
        dataset.read( hydrogenMassFraction, PredType::NATIVE_FLOAT );
        delete group;


        // open the file that gives the HI fraction
        if ( userOptions.additionalOptions.empty() )
            throwError( "You need to give the name of the HDF5 files giving the Urchin HI fraction if you want to compute the HI mass distribution. This is inserted using the '--options' program option." );
        std::string h1FileName = gadgetHeader.filename( userOptions.additionalOptions[0], fileIndex );
        H5File *file2 = new H5File( h1FileName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );
        group = new Group( file2->openGroup( "/PartType0" ) );
        DataSet datasetMolecular = group->openDataSet("MolecularHydrogenMassFraction");
        DataSet datasetHydrogen = group->openDataSet("HydrogenOneFraction");

        float *molecularMassFraction = new float[ gadgetHeader.npart[0] ];
        datasetMolecular.read( molecularMassFraction, PredType::NATIVE_FLOAT );
        float *hydrogenOneFraction   = new float[ gadgetHeader.npart[0] ];
        datasetHydrogen.read( hydrogenOneFraction, PredType::NATIVE_FLOAT );
        delete group;
        delete file2;


        // get the HI mass
        float *scalar = readData->scalar();          // returns a pointer to the particle scalar properties array
        float *weights = readData->weight();         // returns a pointer to the particle weights array
        size_t dataOffset = (*numberParticlesRead);  // the offset in the scalar array from where to start reading the new values
        for (size_t i=0; i<size_t(gadgetHeader.npart[0]); ++i)
        {
            size_t index1 = dataOffset + i;
            size_t index2 = index1 * NO_SCALARS + noScalarsRead;
            scalar[index2] = weights[index1] * hydrogenMassFraction[i] * (float(1.)-molecularMassFraction[i]) * hydrogenOneFraction[i];
        }

        delete[] hydrogenMassFraction;
        delete[] molecularMassFraction;
        delete[] hydrogenOneFraction;
        noScalarsRead += 1;
    }

    delete file;
    for (int i=0; i<6; ++i)     // update the number of read particles
        (*numberParticlesRead) += gadgetHeader.npart[i];
}



/*! Function that counts how many particles are in a Gadget file snapshot that was saved in multiple files. Works only for HDF5 files.
 */
void HDF5_countGadgetParticleNumber(std::string filenameRoot,
                                    int const noFiles,
                                    int const verboseLevel,
                                    size_t numberTotalParticles[])
{
    for (int i=0; i<6; ++i)
        numberTotalParticles[i] = 0;

    // loop over all the files and count the number of particles
    for (int i=0; i<noFiles; ++i)
    {
        Gadget_header gadgetHeader;
        std::string fileName = gadgetHeader.filename( filenameRoot, i );

        // read the header
        HDF5_readGadgetHeader( fileName, &gadgetHeader );

        // add the particle numbers
        for (int j=0; j<6; ++j)
            numberTotalParticles[j] += gadgetHeader.npart[j];
    }

    MESSAGE::Message message( verboseLevel );
    message << "The data is in " << noFiles << " files and contains the following number of particles: " << numberTotalParticles[0] << " + "  << numberTotalParticles[1] << " + "  << numberTotalParticles[2] << " + "  << numberTotalParticles[3] << " + "  << numberTotalParticles[4] << " + "  << numberTotalParticles[5] << " .\n" << MESSAGE::Flush;
}




/*! Read the gadget header for one snapshot and intialize some of the values in the userOptions class.
*/
void HDF5_initializeGadget(std::string filename,
                           Read_data<float> *readData,
                           User_options *userOptions,
                           Gadget_header *gadgetHeader,
                           size_t *noParticles)
{
    MESSAGE::Message message( userOptions->verboseLevel );
    std::string fileName = filename;
    bool singleFile = true;

    // check to see if there is only one input file or several
    if ( not bfs::exists(fileName) ) //if this is true, than the input is in several files
    {
        fileName = gadgetHeader->filename( filename, 0 );
        singleFile = false;
    }


    // now read the actual values of the gadget header
    HDF5_readGadgetHeader( fileName, gadgetHeader );


    // check if to set the box coordinates from the information in the header
    if ( not userOptions->userGivenBoxCoordinates )
    {
        for (size_t i=0; i<NO_DIM; ++i)
        {
            userOptions->boxCoordinates[2*i] = 0.;                    // this is the left extension of the full box
            userOptions->boxCoordinates[2*i+1] = gadgetHeader->BoxSize;// right extension of the full box
        }
    }
    else
        message << "The box coordinates were set by the user using the program options. The program will keep this values and will NOT use the box length information from the Gadget file!" << MESSAGE::Flush;
#ifdef WOJTEK
    if ( userOptions->additionalOptions.size()!=0 ) //if inserted an option
        gadgetHeader->num_files = atoi( userOptions->additionalOptions[0].c_str() );
    gadgetHeader->print();
#endif


    // get the total number of particles to be read from the file/files
    size_t numberTotalParticles[6];
    if ( singleFile )    // if single file
    {
        gadgetHeader->num_files = 1;
        for (int i=0; i<6; ++i)
            numberTotalParticles[i] = gadgetHeader->npart[i];
    }
    else                                // for multiple files read the number of particles in each file and keep track of that
        HDF5_countGadgetParticleNumber( filename, gadgetHeader->num_files, userOptions->verboseLevel, numberTotalParticles );

    // from that keep only the particles specified by the user (not all the particle species may be of interest)
    *noParticles = 0;   //total number of particles to be read from file
    for (int i=0; i<6; ++i)
    {
        if ( not userOptions->readParticleSpecies[i] )
            numberTotalParticles[i] = 0;
        *noParticles += numberTotalParticles[i];
    }
    message << "Reading " << *noParticles << " particle data from the input file. These particles are made from different types of number: ";
    for (int i=0; i<6; ++i)
        message << numberTotalParticles[i] << " (" << gadgetHeader->npart[i] << ") + ";
    message << " .\n" << MESSAGE::Flush;



    // allocate memory for the particle data
    if ( userOptions->readParticleData[0] )
        readData->position( *noParticles );  // particle positions
    if ( userOptions->readParticleData[1] )
        readData->weight( *noParticles );    // particle weights (weights = particle masses from the snapshot file)
#ifdef VELOCITY
    if ( userOptions->readParticleData[2] )
        readData->velocity( *noParticles );  // particle velocities
#endif
    // allocate memory for scalar quantities
#ifdef SCALAR
    int noScalars = 0;
    for (size_t i=3; i<userOptions->readParticleData.size(); ++i)
        if ( userOptions->readParticleData[i] )
            noScalars += 1;
    if ( noScalars>0 )
        readData->scalar( *noParticles );  // particle scalar quantity
#endif



    // check that the 'readParticleData' and 'readParticleSpecies' options make sense
    if ( not userOptions->readParticleData[0] )
        throwError( "The program needs the particle position information to be able to interpolate the fields on a grid. Please add '1' to the integer number giving the data blocks to be read from the input Gadget snapshot." );
    if ( not userOptions->readParticleData[1] )
    {
        MESSAGE::Warning warning( userOptions->verboseLevel );
        warning << "You selected not to read the Gadget particle masses. This means that the program will treat all particles as having the same weight (mass)." << MESSAGE::EndWarning;
    }
    if ( (*noParticles)<=0 )
        throwError( "Please select again the particle species that you would like to read from the Gadget file. There are no particles in the current selection!" );
#ifdef SCALAR
    if ( noScalars>NO_SCALARS )
        throwError( "You asked for too many scalars to be read from the file using the input data option. Either ask for less scalar quantities or increase the number of scalar field components using the Makefile option '-DNO_SCALARS'. You asked for ", noScalars, " scalar components, but 'NO_SCALARS' is only ", NO_SCALARS, "." );
#endif
}



/*! The function that reads the data from single or multiple HDF5 files.
*/
void HDF5_readGadgetFile(std::string filename,
                         Read_data<float> *readData,
                         User_options *userOptions)
{
    MESSAGE::Message message( userOptions->verboseLevel );

    // get the general properties about the data
    Gadget_header gadgetHeader;
    size_t noParticles = 0;
    HDF5_initializeGadget( filename, readData, userOptions, &gadgetHeader, &noParticles );


    // read the data from the files
    size_t numberParticlesRead = 0;
    std::string fileName;
    for (int i=0; i<gadgetHeader.num_files; ++i)
    {
        fileName = gadgetHeader.filename( filename, i );
        message << "Reading GADGET snapshot file '" << fileName << "' which is file " << i+1 << " of " << gadgetHeader.num_files << " files...\n" << MESSAGE::Flush;

        HDF5_readGadgetData( fileName, readData, *userOptions, i, &numberParticlesRead );
    }
}






#endif





