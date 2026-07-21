#ifndef NOAHMP_IO_H
#define NOAHMP_IO_H

#include <NoahArray.H>
#include <vector>

extern "C" {
    struct NoahmpIO_type_fi {

        int *ids, *ide, *jds, *jde, *kds, *kde;
        int *ims, *ime, *jms, *jme, *kms, *kme;
        int *its, *ite, *jts, *jte, *kts, *kte;
        int *xstart, *xend, *ystart, *yend;
        int *nsoil, *nsnow;
        int *itimestep, *ntime;
        int *rank, *blkid, *level;
        int *comm;

        double* XLAT = nullptr;       // latitude [rad]
        double* WSLAKEXY = nullptr;   // lake water storage [mm]
        double* U_PHY = nullptr;
        double* T_PHY = nullptr;
        double* V_PHY = nullptr;
        double* QV_CURR = nullptr;
        double* HFX = nullptr;
        double* LH = nullptr;
        double* SWDOWN = nullptr;
        double* GLW = nullptr;
        double* TSK = nullptr;
        double* EMISS = nullptr;
        double* ALBSFCDIRXY = nullptr;
        double* ALBSFCDIFXY = nullptr;
        double* COSZEN = nullptr;
        double* P8W = nullptr;
        double* TAU_EW = nullptr;
        double* TAU_NS = nullptr;

        NoahmpIO_type_fi(        
          int* ids, int* ide, int* jds, int* jde, int* kds, int* kde,
          int* ims, int* ime, int* jms, int* jme, int* kms, int* kme,
          int* its, int* ite, int* jts, int* jte, int* kts, int* kte,
          int* xstart, int* xend, int* ystart, int* yend, int* itimestep,
          int* ntime, int* nsoil, int* nsnow, int* rank,
          int* blkid, int* level, int *comm
        ):
          ids(ids), ide(ide), jds(jds), jde(jde), kds(kds), kde(kde),
          ims(ims), ime(ime), jms(jms), jme(jme), kms(kms), kme(kme),
          its(its), ite(ite), jts(jts), jte(jte), kts(kts), kte(kte),
          xstart(xstart), xend(xend), ystart(ystart), yend(yend), itimestep(itimestep),
          ntime(ntime), nsoil(nsoil), nsnow(nsnow), rank(rank), blkid(blkid), 
          level(level), comm(comm) {}
    };
}

class NoahmpIO_type {

    public:
        int ids, ide, jds, jde, kds, kde;
        int ims, ime, jms, jme, kms, kme;
        int its, ite, jts, jte, kts, kte;
        int xstart, xend, ystart, yend;
        int nsoil, nsnow;
        int itimestep, ntime;
        int rank, blkid, level;
        int comm;

        NoahArray2D<double> XLAT;
        NoahArray2D<double> WSLAKEXY;
        NoahArray3D<double> U_PHY;
        NoahArray3D<double> T_PHY;
        NoahArray3D<double> V_PHY;
        NoahArray3D<double> QV_CURR;
        NoahArray2D<double> HFX;
        NoahArray2D<double> LH;
        NoahArray2D<double> SWDOWN;
        NoahArray2D<double> GLW;
        NoahArray2D<double> TSK;
        NoahArray2D<double> EMISS;
        NoahArray3D<double> ALBSFCDIRXY;
        NoahArray3D<double> ALBSFCDIFXY;
        NoahArray2D<double> COSZEN;
        NoahArray3D<double> P8W;
        NoahArray2D<double> TAU_EW, TAU_NS;

        NoahmpIO_type():
          fptr(&ids, &ide, &jds, &jde, &kds, &kde,
               &ims, &ime, &jms, &jme, &kms, &kme,
               &its, &ite, &jts, &jte, &kts, &kte,
               &xstart, &xend, &ystart, &yend, &itimestep,
               &ntime, &nsoil, &nsnow, &rank, &blkid, &level, &comm) {}

        void ScalarInitDefault();
        void VarInitDefault();
        void InitMain();
        void ReadNamelist();
        void ReadTable();
        void ReadLandHeader();
        void ReadLandMain();
        void DriverMain();
        void WriteLand(int filenum);

    private:
        NoahmpIO_type_fi fptr;
};

class NoahmpIO_vector : public std::vector<NoahmpIO_type> {
public:
    void resize(size_t size, size_t level);

};

#endif // NOAHMP_IO_H
