// #include "/group/nps/hhuang/analysis/MyHeader/Analysis.h"
#include "/group/nps/hhuang/analysis/MyHeader/MyDVCSDB.h"
#include "/group/nps/hhuang/software/NPS_SOFT/TDVCSDB.h"
#include <filesystem>
using namespace std;
namespace fs = std::filesystem;
ofstream logfile;
Int_t GetNumberOfCycles(int kinc_param1, const string& kinc_param2, int target_flag);

void checkStatistics(int kinc_param1, string kinc_param2, int target_flag = -1)
{
    gSystem->Load("/group/nps/hhuang/software/NPS_SOFT/libDVCS.so");
    TDVCSDB *db = new TDVCSDB("dvcs", "clrlpc", 3306, "hhuang", "");

    // Get the run range of each kinematics____________________________________
    vector<int> minrun, maxrun;
    GetRunRange(kinc_param1, kinc_param2, minrun, maxrun);
    TString kinName = Form("x%d_%s", kinc_param1, kinc_param2.c_str());

    // using the run list of cycles to count total number of cycles
    Int_t nSubKine = GetNumberOfCycles(kinc_param1, kinc_param2, target_flag);

    // Files output to check pi0 statistics
    TFile *outfile = new TFile(Form("output/%s_%d_wf.root", kinName.Data(), target_flag), "recreate");
    ofstream fsummary;
    fsummary.open(Form("output/%s_%d_wf.txt", kinName.Data(), target_flag));

    // Files output to check background ratio
    ofstream fsummary_bkg;
    fsummary_bkg.open(Form("output/%s_%d_wf_bkg.txt", kinName.Data(), target_flag));
    
    // Histograms to check statistics
    TH1F *h_pi0Stat[nSubKine][11];
    for(int isk = 0; isk < nSubKine; isk++){
        Double_t x1 = minrun[0]-0.5;
        Double_t x2 = maxrun[maxrun.size()-1]+0.5;
        Int_t nbins = x2-x1;
        for(int i = 0; i < 11; i++){
            h_pi0Stat[isk][i] = new TH1F(Form("h_pi0Stat_cycle%d_%d", isk, i), "#pi^{0} statistics;Run number;Number of #pi^{0}", nbins, x1, x2);
        }
    }

    // Get the pi0 histograms from rootfiles, fit and check the statistis
    ifstream fRunList;
    vector<int> runList;
    vector<int> nSegList;
    string line;

    fsummary<<"This is a list of pi0 statistics with different energy cut on photons"<<endl;
    fsummary<<"Run#  0.5GeV  0.6GeV  0.7GeV  0.8GeV  0.9GeV  1.0GeV  1.1GeV  1.2GeV  1.3GeV  1.4GeV  1.5GeV"<<endl;
    fsummary<<endl;

    fsummary_bkg<<"This is a list of background ratio with different energy cut on photons"<<endl;
    fsummary_bkg<<"Run#  0.5GeV  0.6GeV  0.7GeV  0.8GeV  0.9GeV  1.0GeV  1.1GeV  1.2GeV  1.3GeV  1.4GeV  1.5GeV"<<endl;
    fsummary_bkg<<endl;

    for(int isk = 0; isk < nSubKine; isk++){ // loop of cycles
        // Get the list of runs
        fRunList.open(Form("/group/nps/hhuang/analysis/DVCS_NPS2023/DVCS_analysis/pi0Calib_wf/calibList/ListPerCycle/%s_%d_cycle%d.txt", kinName.Data(), target_flag, isk));
        if(!fRunList){
            cout<<"ERROR: Can't find the list of runs for calibration!!!"<<endl;
            return;
        }

        while (getline(fRunList, line)){
            istringstream iss(line);
            Int_t temp_run;
            Int_t temp_nSeg;
            if(iss >> temp_run >> temp_nSeg){
                runList.push_back(temp_run);
                nSegList.push_back(temp_nSeg);
            }
        }
        fRunList.close(); // close the text file of run list

        Int_t nRun = runList.size();
        Int_t first_run = runList[0];
        Int_t last_run = runList[nRun-1];

        // cout<<"This is cycle "<<i<<endl;
        // for(int irun = 0; irun < nRun; irun++){
        //     cout<<"Run: "<<runList[irun]<<"; nSeg: "<<nSegList[irun]<<endl;
        // }

        // Get histograms from each rootfile, fit and check the statistics
        TFile *infile;
        TString infileDir = "/group/nps/hhuang/analysis/DVCS_NPS2023/DVCS_analysis/pi0Calib_wf/calibTree";
        TString filename_prefix = "prodTree_wf";
        TString filename;
        TH1F *h_pi0M[11] = {nullptr};
        Int_t nEntries[11]; // To check the histogram after merging
        Double_t sumNpi0[11]; // To sum over runs of the cycle for signal
        for(int i = 0; i < 11; i++) sumNpi0[i] = 0;
        Double_t sumNbkg[11]; // To sum over runs of the cycle for background
        for(int i = 0; i < 11; i++) sumNbkg[i] = 0;

        fsummary<<"++++++++++ This is cycle "<<isk<<" ++++++++++"<<endl;
        fsummary_bkg<<"++++++++++ This is cycle "<<isk<<" ++++++++++"<<endl;
        for(int irun = 0; irun < nRun; irun++){ // loop of runs in that cycle
            if(runList[irun]==4745||runList[irun]==4803||runList[irun]==4713||runList[irun]==4715) continue;
            for(int i = 0; i < 11; i++) nEntries[i] = 0;

            for(int iseg = 0; iseg < nSegList[irun]; iseg++){ // loop the segment of each run
                filename = Form("%s/%s/%s_%d_%d.root", infileDir.Data(), kinName.Data(), filename_prefix.Data(), runList[irun], iseg);
                infile = TFile::Open(filename.Data());
                if (!infile || infile->IsZombie()) {
                    cerr << "Error: cannot open file for Run "<<runList[irun]<<", segment "<<iseg<<endl;
                    delete infile;
                    infile = nullptr;
                    continue;
                }

                // Get and merge histograms
                for(int i = 0; i < 11; i++){
                    TH1F *h_temp = (TH1F*)infile->Get(Form("h_pi0M_%d", i)); // Get the histogram of pi0 mass
                    if (!h_temp){
                        cerr << "Warning: histogram h_pi0M_" << i << " not found in " << filename << endl;
                        continue;
                    }
                    nEntries[i]+=h_temp->GetEntries(); // to check if the histograms are merged correctly
                    if(!h_pi0M[i]){ // Clone if it is from the first segment
                        h_pi0M[i] = (TH1F*)h_temp->Clone(Form("h_pi0M_run%d_%d", runList[irun], i));
                        h_pi0M[i]->SetDirectory(0);
                    }
                    else h_pi0M[i]->Add(h_temp, 1);
                }

                delete infile;
                infile = nullptr;
            } // loop of segments

            for(int i = 0; i < 11; i++){ // check if the histograms are merged correctly
                if(h_pi0M[i]->GetEntries() != nEntries[i]) cout<<"Error: histogram "<<i<<" of Run "<<runList[irun]<<" was not added correctly"<<endl;
            }

            fsummary<<runList[irun]<<" ";
            fsummary_bkg<<runList[irun]<<" ";
            for(int i = 0; i < 11; i++){ // Fit to get the statistics
                Double_t meanfit = h_pi0M[i]->GetBinCenter(h_pi0M[i]->GetMaximumBin());
                Double_t minfit = meanfit - 0.04;
                Double_t maxfit = meanfit + 0.04;
                TF1 *f_fit = new TF1(Form("f_fit_cycle%d_%d", isk, i), "gausn(0)+pol1(3)", minfit, maxfit);
                f_fit->SetParLimits(1, 0.1, 0.15);
                f_fit->SetParLimits(2, 0.00001, 0.03);
                h_pi0M[i]->Fit(f_fit, "R");

                Float_t N = f_fit->GetParameter(0)/h_pi0M[i]->GetBinWidth(1);
                Float_t N_err = f_fit->GetParError(0)/h_pi0M[i]->GetBinWidth(1);
                Float_t mean = f_fit->GetParameter(1);
                Float_t sigm = f_fit->GetParameter(2);

                Int_t ibin = h_pi0Stat[isk][i]->FindBin(runList[irun]);
                h_pi0Stat[isk][i]->SetBinContent(ibin, N);
                h_pi0Stat[isk][i]->SetBinError(ibin, N_err);

                sumNpi0[i]+=N;
                fsummary<<N<<" ";

                // calculate background ratio
                TF1 *f_bkg = new TF1(Form("f_bkg_cycle%d_%d", isk, i), "pol1(0)", minfit, maxfit);
                f_bkg->SetParameters(f_fit->GetParameter(3), f_fit->GetParameter(4));
                Double_t bkg_integral = f_bkg->Integral(mean-3*sigm, mean+3*sigm)/h_pi0M[i]->GetBinWidth(1);
                sumNbkg[i]+=bkg_integral;
                fsummary_bkg<<bkg_integral/N<<" "; // background ratio

                outfile->cd();
                h_pi0M[i]->Write();
                delete f_fit;
                delete f_bkg;
            }
            fsummary<<endl;
            fsummary_bkg<<endl;

            for(int i = 0; i < 11; i++){
                delete h_pi0M[i];
                h_pi0M[i] = nullptr;
            }
        } // loop of runs
        fsummary<<endl;
        fsummary_bkg<<endl;
        fsummary<<"Total ";
        fsummary_bkg<<"Total ";
        for(int i = 0; i < 11; i++) fsummary<<sumNpi0[i]<<" ";
        for(int i = 0; i < 11; i++) fsummary_bkg<<sumNbkg[i]/sumNpi0[i]<<" ";
        fsummary<<endl;
        fsummary<<endl;
        fsummary_bkg<<endl;
        fsummary_bkg<<endl;
        runList.clear();
        nSegList.clear();
    } // loop of cycles
    fsummary.close();
    fsummary_bkg.close();

    outfile->cd();
    for(int isk = 0; isk < nSubKine; isk++){
        for(int i = 0; i < 11; i++){
            h_pi0Stat[isk][i]->Write();
        }
    }
}

Int_t GetNumberOfCycles(int kinc_param1, const string& kinc_param2, int target_flag)
{   
    const string dir = "/group/nps/hhuang/analysis/DVCS_NPS2023/DVCS_analysis/pi0Calib_wf/calibList/ListPerCycle/";
    string pattern_str = "x"+to_string(kinc_param1)+"_"+kinc_param2+"_"+to_string(target_flag)+"_cycle(\\d+)\\.txt";
    regex pattern(pattern_str);

    int max_cycle = -1;

    for (const auto& entry : fs::directory_iterator(dir)) {
        if (!entry.is_regular_file()) continue;

        string fname = entry.path().filename().string();
        smatch match;

        if(regex_match(fname, match, pattern)) {
            int cycle = stoi(match[1]);
            max_cycle = max(max_cycle, cycle);
        }
    }

    return (max_cycle >= 0) ? max_cycle + 1 : 0;
}

// void MakeListPerCycle()
// {
//     // 54 Kinematics
//     const Int_t nKine = 54;
//     Int_t array_par1[nKine] = {
//         // KinC_x25 (5)
//         25, 25, 25, 25, 25,
//         // KinC_x36 (16)
//         36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36,
//         // KinC_x50 (19)
//         50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50,
//         // KinC_x60 (14)
//         60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60
//     };
//     string array_par2[nKine] = {
//         // KinC_x25
//         "1_1", "1_2", "3", "4_1", "4_2",
//         // KinC_x36
//         "1", "2_1", "2_2", "2p", "2pp", "3_1", "3_2", "4", "5_1", "5_2", "5_3", "5_4", "5p", "6_1", "6_2", "6_3",
//         // KinC_x50
//         "0_1", "0_2", "0a_1", "0a_2", "0b_1", "0b_2", "1_1", "1_2", "1p", "2_1", "2_2", "2_3", "2p", "2pp", "3_1", "3_2", "3_3", "3p", "3pp",
//         // KinC_x60
//         "1", "2b", "2_1", "2_2", "2p", "3_1", "3_2", "3_3", "3p", "3a", "3b", "4a_1", "4a_2", "4b"
//     };

//     system("rm -f ListPerCycle/log_cycleList.txt");
//     logfile.open("ListPerCycle/log_cycleList.txt");
    
//     for(int i=0; i<54; i++){
//         MakeList(array_par1[i], array_par2[i]);
//     }
//     cout<<"Finished making all calibration lists"<<endl;
// }