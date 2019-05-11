/***********************************************************************
/
/  GRID CLASS (IDENTIFY THE MULTI-SPECIES FIELDS)
/
/  written by: Greg Bryan
/  date:       October, 1996
/  modified1:
/
/  PURPOSE:
/
/  NOTE:
/
************************************************************************/
 
#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "fortran.def"
 
/* function prototypes */
 
int FindField(int f, int farray[], int n);
 
 
 
int grid::IdentifySpeciesFields(int &DeNum, int &HINum, int &HIINum,
				int &HeINum, int &HeIINum, int &HeIIINum,
                                int &HMNum, int &H2INum, int &H2IINum,
                                int &DINum, int &DIINum, int &HDINum)
{
 
  DeNum = HINum = HIINum = HeINum = HeIINum = HeIIINum = 0;
  HMNum = H2INum = H2IINum = DINum = DIINum = HDINum = 0;
 
  /* Find Fields for the 6-species model. */
 
  DeNum    = FindField(ElectronDensity, FieldType, NumberOfBaryonFields);
  HINum    = FindField(HIDensity, FieldType, NumberOfBaryonFields);
  HIINum   = FindField(HIIDensity, FieldType, NumberOfBaryonFields);
  HeINum   = FindField(HeIDensity, FieldType, NumberOfBaryonFields);
  HeIINum  = FindField(HeIIDensity, FieldType, NumberOfBaryonFields);
  HeIIINum = FindField(HeIIIDensity, FieldType, NumberOfBaryonFields);
 
  /* Error if any not found. */
 
  if (DeNum < 0 || HINum < 0 || HIINum < 0 || HeINum < 0 || HeIINum < 0 ||
      HeIIINum < 0) {
    ENZO_VFAIL("De=%"ISYM", HI=%"ISYM", HII=%"ISYM", HeI=%"ISYM", HeII=%"ISYM", HeIII=%"ISYM"\n",
	    DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum)
  }
 
  /* Find Fields for the 9-species model. */
 
  if (MultiSpecies > 1) {
    HMNum    = FindField(HMDensity, FieldType, NumberOfBaryonFields);
    H2INum   = FindField(H2IDensity, FieldType, NumberOfBaryonFields);
    H2IINum  = FindField(H2IIDensity, FieldType, NumberOfBaryonFields);
 
    if (HMNum < 0 || H2INum < 0 || H2IINum < 0) {
      ENZO_FAIL("H2 related field missing.\n");
    }
 
  }
 
  /* Find Fields for the 12-species model. */
 
  if (MultiSpecies > 2) {
    DINum   = FindField(DIDensity, FieldType, NumberOfBaryonFields);
    DIINum  = FindField(DIIDensity, FieldType, NumberOfBaryonFields);
    HDINum  = FindField(HDIDensity, FieldType, NumberOfBaryonFields);
 
    if (DINum < 0 || DIINum < 0 || HDINum < 0) {
      ENZO_FAIL("HD related field missing.\n");

    }
 
  }
 
  return SUCCESS;
}

 int grid::IdentifySpeciesFieldsChem(int &WaterNum, int &ONum, int &OHNum, 
                                     int &O2Num, int &OplusNum,
                                     int &OHplusNum, int &H2OplusNum, 
                                     int &H3OplusNum, int &O2plusNum, 
                                     int &CplusNum, int &CNum, 
                                     int &CHNum, int &CH2Num, 
                                     int &CH3Num, int &CH4Num, 
                                     int &CONum, int &COplusNum, int &CO2Num,
                                     int &SNum, int &SplusNum, int &SiNum, 
                                     int &SiplusNum, int &CHplusNum, int &CH2plusNum, 
                                     int &SiONum, int &H3plusNum)
{
   if (withWater){

       WaterNum = ONum = OHNum = O2Num = OplusNum = OHplusNum = H2OplusNum = H3OplusNum = O2plusNum = 0;
       CplusNum = CNum = CHNum = CH2Num = CH3Num = CH4Num = CONum = COplusNum = CO2Num = 0;
       SNum = SplusNum = SiNum = SiplusNum = CHplusNum = CH2plusNum = SiONum = H3plusNum = 0; 

       WaterNum = FindField(WaterDensity, FieldType, NumberOfBaryonFields);
       ONum = FindField(ODensity, FieldType, NumberOfBaryonFields);
       OHNum = FindField(OHDensity, FieldType, NumberOfBaryonFields);
       O2Num = FindField(O2Density, FieldType, NumberOfBaryonFields);
       OplusNum = FindField(OplusDensity, FieldType, NumberOfBaryonFields);
       OHplusNum = FindField(OHplusDensity, FieldType, NumberOfBaryonFields);
       H2OplusNum = FindField(H2OplusDensity, FieldType, NumberOfBaryonFields);
       H3OplusNum = FindField(H3OplusDensity, FieldType, NumberOfBaryonFields);
       O2plusNum = FindField(O2plusDensity, FieldType, NumberOfBaryonFields);
       CplusNum = FindField(CplusDensity, FieldType, NumberOfBaryonFields);
       CNum = FindField(CDensity, FieldType, NumberOfBaryonFields);
       CHNum = FindField(CHDensity, FieldType, NumberOfBaryonFields);
       CH2Num = FindField(CH2Density, FieldType, NumberOfBaryonFields);
       CH3Num = FindField(CH3Density, FieldType, NumberOfBaryonFields);
       CH4Num = FindField(CH4Density, FieldType, NumberOfBaryonFields);
       CONum = FindField(CODensity, FieldType, NumberOfBaryonFields);
       COplusNum = FindField(COplusDensity, FieldType, NumberOfBaryonFields);
       CO2Num = FindField(CO2Density, FieldType, NumberOfBaryonFields);

       if (WaterNum < 0 || ONum < 0 || OHNum < 0 || O2Num < 0 || OplusNum < 0 ||
           OHplusNum < 0 || H2OplusNum < 0 || H3OplusNum < 0 || O2plusNum < 0 || CplusNum < 0 || CNum < 0 ||
           CHNum < 0 || CH2Num < 0 || CH3Num < 0 || CH4Num < 0 ||
           CONum < 0 || COplusNum < 0 || CO2Num < 0)
          {
               ENZO_VFAIL("Water=%"ISYM", O=%"ISYM", OH=%"ISYM", O2=%"ISYM", Oplus=%"ISYM", OHplus=%"ISYM", H2Oplus=%"ISYM", "
              "H3Oplus=%"ISYM", O2plus=%"ISYM", Cplus=%"ISYM", C=%"ISYM", CH=%"ISYM", CH2=%"ISYM", CH3=%"ISYM", CH4=%"ISYM", "
              "CO=%"ISYM", COplus=%"ISYM", CO2=%"ISYM"\n", WaterNum, ONum, OHNum, O2Num, OplusNum, OHplusNum, H2OplusNum,
              H3OplusNum, O2plusNum, CplusNum, CNum, CHNum, CH2Num, CH3Num, CH4Num, CONum, COplusNum, CO2Num);
           }

       if (water_rates == 3) {
         SNum       = FindField(SDensity, FieldType, NumberOfBaryonFields);
         SplusNum   = FindField(SplusDensity, FieldType, NumberOfBaryonFields);
         SiNum      = FindField(SiDensity, FieldType, NumberOfBaryonFields);
         SiplusNum  = FindField(SiplusDensity, FieldType, NumberOfBaryonFields);
         CHplusNum  = FindField(CHplusDensity, FieldType, NumberOfBaryonFields);
         CH2plusNum = FindField(CH2plusDensity, FieldType, NumberOfBaryonFields);
         SiONum     = FindField(SiODensity, FieldType, NumberOfBaryonFields);
         H3plusNum  = FindField(H3plusDensity, FieldType, NumberOfBaryonFields);

         if (SNum < 0 || SplusNum < 0 || SiNum < 0 || SiplusNum < 0 || CHplusNum < 0 || CH2plusNum < 0 || 
             SiONum < 0 ||  H3plusNum < 0) {
           ENZO_FAIL("Bialy Water related field missing.\n");
         }
      }
    }
    return SUCCESS;
}

