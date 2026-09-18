#include "TSystem.h"
#include "TMatrix.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TFile.h"
#include "TCutG.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TGaxis.h"
#include "TGraph.h"
#include "TLorentzVector.h"
#include "TNtuple.h"
#include "TLegend.h"
#include "TLine.h"
#include "TRandom3.h"
#include "TStyle.h"
#include <iostream>
#include <sstream>
// #include "TDVCSGlobal.h"
#include "TTreeIndex.h"
#include "TChainIndex.h"
#include "TCaloEvent.h"
#include "TCaloGeometry.h"
#include "TCaloBase.h"
#include "TDVCSEvent.h"

using namespace std;

Double_t m_p = 0.938272013; // proton mass [GeV]
Double_t m_pi0 = 0.1349766; // pi0 mass [GeV]

const Int_t ncol = 30;          // number of columns
const Int_t nrow = 36;          // number of rows
const Int_t nblk = ncol*nrow; // number of blocks of NPS

Int_t GetColumn(int iblk) // iblk should be in current NPS numbering scheme
{
    Int_t icol = iblk%ncol;
    return icol;
}

Int_t GetRow(int iblk) // iblk should be in current NPS numbering scheme
{
    Int_t icol = iblk%ncol;
    Int_t irow = (iblk-icol)/ncol;
    return irow;
}


// Functions to convert the numbering scheme of NPS
Int_t bnConv_OldToNew(int ibn_old) // Convert the PMT number to the current version
{
  const Int_t ncol = 30; // number of columns
  const Int_t nrow = 36; // number of rows
  Int_t irow = ibn_old % nrow;
  Int_t icol = 29 - ((ibn_old - irow) / nrow);
  Int_t ibn_new = ncol * irow + icol;

  return ibn_new;
}

Int_t bnConv_NewToOld(int ibn_new) // Convert the PMT number to the simulation version
{
  const Int_t ncol = 30; // number of columns
  const Int_t nrow = 36; // number of rows
  Int_t icol = 29 - (ibn_new % ncol);
  Int_t irow = (ibn_new - ibn_new % ncol) / ncol;
  Int_t ibn_old = nrow * icol + irow;

  return ibn_old;
}
/*
Int_t bnConv_FidToFull(int ibn_fid) // only used to convert the fiducial to simulation numbering scheme
{
  // Fiducial row and column number
  Int_t nrow_fid = 32;
  Int_t ncol_fid = 22;
  Int_t irow_fid = ibn_fid%nrow_fid;
  Int_t icol_fid = (ibn_fid-irow_fid)/nrow_fid;

  Int_t irow_sim = irow_fid+2;
  Int_t icol_sim = icol_fid+2;
  Int_t ibn_sim = icol_sim*nrow+irow_sim;

  return ibn_sim;
}
*/