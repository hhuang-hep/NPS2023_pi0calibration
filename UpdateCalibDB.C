#include "/group/nps/hhuang/software/NPS_SOFT/TDVCSDB.h"
#include "TSystem.h"
#include "TSystemDirectory.h"
#include "TList.h"
#include <string>
#include <fstream>
#include <iostream>
using namespace std;

bool findCoefFile(int run, TString& foundFile)
{
    TString fname = Form("coef_pi0Calib_%d.txt", run);

    TSystemDirectory base("Result", "./Result");
    TList* list = base.GetListOfFiles();
    if (!list) return false;

    TIter next(list);
    TObject* obj;

    while ((obj = next())) {
        TString d = obj->GetName();
        if (d == "." || d == "..") continue;

        TString fullPath = "./Result/" + d + "/" + fname;
        if (!gSystem->AccessPathName(fullPath)) {
            foundFile = fullPath;
            return true;
        }
    }
    return false;
}

void UpdateCalibDB(int kinc_param1, string kinc_param2){ // update coefficients with NPS numbering scheme

    // Connect to the database
    gSystem->Load("/group/nps/hhuang/software/NPS_SOFT/libDVCS.so");
    TDVCSDB *db=new TDVCSDB("dvcs","clrlpc",3306,"hhuang","");

    // Tables to change
    char const *T_Pi0Coef = "CALO_calib_Pi0Coef";

    string kinc_param1_str = to_string(kinc_param1);
    string kinName = "x"+kinc_param1_str+"_"+kinc_param2;
    TString firstValidFile;
    TString lastValidFile;
    vector<int> skippedRuns;
    
    for(int irun = 1572; irun <= 7012; irun++){
        const char** val = db->GetEntry_s("Run_param_KineNameSimu", irun); // simulation kinematics name
        string str(*val);
        if(str == "KinC_"+kinName){
            TString coefFile; // coefficient file name
            bool found = findCoefFile(irun, coefFile); // search for the coefficient file in the Result folder

            if (found) { // update the last valid coefficient file when the file of this run is found
                lastValidFile = coefFile;
                cout<<"[run "<<irun<<"] use NEW file:"<<coefFile<<endl;
                if(firstValidFile.IsNull()){
                    firstValidFile = coefFile;
                }
            }
            else { // using the last valid coefficient file if the file this run is not found
                if(lastValidFile.IsNull()) { // still before first valid run
                    cout<<"[run "<<irun<<"] No coefficient file found yet!"<<endl;
                    skippedRuns.push_back(irun);
                    continue;
                }
                coefFile = lastValidFile;
                cout<<"[run "<<irun<<"] fallback to: "<<coefFile<<endl;
            }

            ifstream fcoef(coefFile);
            if (!fcoef){
                cerr<<"Failed to open "<<coefFile<<endl;
                return;
            }
            
            Double_t coef[1080];
            for(int iblk = 0; iblk < 1080; iblk++) fcoef>>coef[iblk];
            // for(int iblk = 0; iblk < 1080; iblk++) cout<<coef[iblk]<<endl;
            db->AddEntry_d(T_Pi0Coef, irun, irun, coef, "pi0 calibration coefficients");
        }
    }

    // In case no coefficient file is found at all
    if(firstValidFile.IsNull()){
        cerr << "ERROR: no valid coefficient file found for "<< kinName << endl;
        return;
    }

    // Handle skipped runs at the beginning
    if(!skippedRuns.empty()){
        ifstream fcoef(firstValidFile);
        if (!fcoef){
            cerr<<"Failed to open "<<firstValidFile<<endl;
            return;
        }

        Double_t coef[1080];
        for(int iblk = 0; iblk < 1080; iblk++) fcoef>>coef[iblk];
        // for(int iblk = 0; iblk < 1080; iblk++) cout<<coef[iblk]<<endl;
        
        for(int irun = 0; irun < skippedRuns.size(); irun++){
            int runNum = skippedRuns[irun];
            cout<<"[run "<<runNum<<"] use FIRST valid file: "<<firstValidFile<<endl;
            db->AddEntry_d(T_Pi0Coef, runNum, runNum, coef, "pi0 calibration coefficients");
        }
    }
}
