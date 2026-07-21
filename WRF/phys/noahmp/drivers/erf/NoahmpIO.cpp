#include <NoahmpIO.H>
#include <NoahArray.H>

extern "C" {
    void NoahmpIOScalarInitDefault_fi(NoahmpIO_type_fi* noahmpio);
    void NoahmpIOVarInitDefault_fi(NoahmpIO_type_fi* noahmpio);
    void NoahmpInitMain_fi(NoahmpIO_type_fi* noahmpio);
    void NoahmpReadNamelist_fi(NoahmpIO_type_fi* noahmpio);
    void NoahmpReadTable_fi(NoahmpIO_type_fi* noahmpio);
    void NoahmpReadLandHeader_fi(NoahmpIO_type_fi* noahmpio);
    void NoahmpReadLandMain_fi(NoahmpIO_type_fi* noahmpio);
    void NoahmpDriverMain_fi(NoahmpIO_type_fi* noahmpio);
    void NoahmpWriteLand_fi(NoahmpIO_type_fi* noahmpio, int* filenum);
    void NoahmpIOTypeVectInit_fi(int* level, int* NBlocks);
}

void NoahmpIO_type::ScalarInitDefault() {
     NoahmpIOScalarInitDefault_fi(&fptr);
};


void NoahmpIO_type::InitMain() {
     NoahmpInitMain_fi(&fptr);
};

void NoahmpIO_type::ReadNamelist() {
     NoahmpReadNamelist_fi(&fptr);
};

void NoahmpIO_type::ReadTable() {
     NoahmpReadTable_fi(&fptr);
};

void NoahmpIO_type::ReadLandHeader() {
     NoahmpReadLandHeader_fi(&fptr);
};

void NoahmpIO_type::ReadLandMain() {
     NoahmpReadLandMain_fi(&fptr);
};

void NoahmpIO_type::DriverMain() {
     NoahmpDriverMain_fi(&fptr);
};

void NoahmpIO_type::WriteLand(int filenum) {
     NoahmpWriteLand_fi(&fptr, &filenum);
};


void NoahmpIO_type::VarInitDefault() {

      NoahmpIOVarInitDefault_fi(&fptr);

      XLAT     = NoahArray2D<double>(fptr.XLAT,     {xstart,ystart}, {xend,yend});
      WSLAKEXY = NoahArray2D<double>(fptr.WSLAKEXY, {xstart,ystart}, {xend,yend});

      T_PHY   = NoahArray3D<double>(fptr.T_PHY,   {xstart,kms,ystart}, {xend,kme,yend});
      U_PHY   = NoahArray3D<double>(fptr.U_PHY,   {xstart,kms,ystart}, {xend,kme,yend});
      V_PHY   = NoahArray3D<double>(fptr.V_PHY,   {xstart,kms,ystart}, {xend,kme,yend});
      QV_CURR = NoahArray3D<double>(fptr.QV_CURR, {xstart,kms,ystart}, {xend,kme,yend});

      HFX = NoahArray2D<double>(fptr.HFX, {xstart,ystart}, {xend,yend});
      LH = NoahArray2D<double>(fptr.LH, {xstart,ystart}, {xend,yend});

      SWDOWN = NoahArray2D<double>(fptr.SWDOWN, {xstart,ystart}, {xend,yend});
      GLW = NoahArray2D<double>(fptr.GLW, {xstart,ystart}, {xend,yend});
      TSK = NoahArray2D<double>(fptr.TSK, {xstart,ystart}, {xend,yend});
      EMISS = NoahArray2D<double>(fptr.EMISS, {xstart,ystart}, {xend,yend});

      ALBSFCDIRXY = NoahArray3D<double>(fptr.ALBSFCDIRXY, {xstart,1,ystart}, {xend,2,yend});
      ALBSFCDIFXY = NoahArray3D<double>(fptr.ALBSFCDIFXY, {xstart,1,ystart}, {xend,2,yend});

      COSZEN = NoahArray2D<double>(fptr.COSZEN, {xstart,ystart}, {xend,yend});
      P8W = NoahArray3D<double>(fptr.P8W, {xstart,kms,ystart}, {xend,kme,yend});

      TAU_EW = NoahArray2D<double>(fptr.TAU_EW, {xstart,ystart}, {xend,yend});
      TAU_NS = NoahArray2D<double>(fptr.TAU_NS, {xstart,ystart}, {xend,yend});
};


void NoahmpIO_vector::resize(size_t size, size_t level) {
     std::vector<NoahmpIO_type>::resize(size);
     int _size = size;
     int _level = level;
     NoahmpIOTypeVectInit_fi(&_level, &_size);
};
