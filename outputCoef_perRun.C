#include "/group/nps/hhuang/analysis/MyHeader/Analysis.h"
#include "/group/nps/hhuang/analysis/MyHeader/MyDB.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include "TVectorD.h"
using namespace std;

// Number of columns used for calibration
const Int_t ncol_inUse = 23;
const Int_t nblk_inUse = ncol_inUse*nrow;

// Check if the block is in the range for calibration (with sim. numbering scheme)
bool passAccCut(int sim_blk_number);

void outputCoef_perRun(TString Kine, int TargetFlag)
{
    const Int_t nIter = 8;

    TString Tar;
    if(TargetFlag == 0) Tar = "LH2";
    else if(TargetFlag == 1) Tar = "LD2";
    else if(TargetFlag == -1) Tar = "LH2_LD2";
    else{
        cout<<"ERROR: Unknown target!!!"<<endl;
        return;
    }

    TString filename = Form("%s_%s", Kine.Data(), Tar.Data());

    // Get the list of runs for calibration
    ifstream fRunList(Form("/group/nps/hhuang/analysis/DVCS_NPS2023/DVCS_analysis/MakeProdList/%s_%d_ProdList.txt", Kine.Data(), TargetFlag));
    vector<int> runList;
    vector<int> nSegList;
    string line;

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

    Int_t nRun = runList.size();
    Int_t first_run = runList[0];
    Int_t last_run = runList[nRun-1];
    fRunList.close();

    // Tree and branches in the rootfiles
    TChain *chain = new TChain("t_prod");;
    Double_t hms_Vx;
    Double_t hms_Vy;
    Double_t hms_Vz;
    TCaloEvent *caloev = new TCaloEvent();

    // Initial settings
    TDVCSEvent *ev = new TDVCSEvent();
    Bool_t hasBlock[1080];

    // Connect to the database to get kinematics variables________________________________________________
    gSystem->Load("/group/nps/hhuang/software/NPS_SOFT/libDVCS.so");
    TDVCSDB *db = new TDVCSDB("dvcs", "clrlpc", 3306, "hhuang", "");

    // Input the coefficients after the last iteration (in NPS numbering scheme)
    ifstream fpi0_coef(Form("/group/nps/hhuang/analysis/DVCS_NPS2023/DVCS_analysis/pi0Calib_wf/Result/%s_pass2_v3/coef_pi0Calib_temp.txt", filename.Data()));
    if(!fpi0_coef || !fpi0_coef.is_open()){
        cout<<"ERROR: Can't find the pi0 coefficients!!!"<<endl;
        return;
    }
    Double_t pi0_coef[1080];
    for(int i = 0; i < 1080; i++) fpi0_coef>>pi0_coef[i]; // in NPS numbering scheme
    
    TFile *outfile = new TFile(Form("Result/%s_pass2_v3/UpdateCoef.root", filename.Data()), "recreate");

    TH1F *h_pi0M = new TH1F("h_pi0M", "M_{#gamma#gamma} (E_{#gamma1} > 1.4 GeV && E_{#gamma2} > 1.4 GeV);M_{#gamma#gamma} [GeV];Counts", 75, 0.05, 0.2);
    h_pi0M->Sumw2();
    TF1 *f_fit = new TF1("f_fit", "gausn(0)+pol1(3)", 0.08, 0.2);
    Float_t N, N_err, mean, mean_err, sigm, sigm_err;
    TH1F *h_pi0MeanvsRun_old = new TH1F("h_pi0MeanvsRun_old", "M_{#pi0} vs Run Number;Run Number;M_{#pi0} [MeV/c^{2}]", last_run-first_run+1, first_run-0.5, last_run+0.5);
    TH1F *h_pi0SigmvsRun_old = new TH1F("h_pi0SigmvsRun_old", "#sigma_{#pi0} vs Run Number;Run Number;#sigma_{#pi0} [MeV/c^{2}]", last_run-first_run+1, first_run-0.5, last_run+0.5);
    TH1F *h_pi0MeanvsRun_new = new TH1F("h_pi0MeanvsRun_new", "M_{#pi0} vs Run Number;Run Number;M_{#pi0} [MeV/c^{2}]", last_run-first_run+1, first_run-0.5, last_run+0.5);
    TH1F *h_pi0SigmvsRun_new = new TH1F("h_pi0SigmvsRun_new", "#sigma_{#pi0} vs Run Number;Run Number;#sigma_{#pi0} [MeV/c^{2}]", last_run-first_run+1, first_run-0.5, last_run+0.5);

    Int_t nevt;
    for(int irun = 0; irun < nRun; irun++){

        for(int iseg = 0; iseg < nSegList[irun]; iseg++){
            chain->Add(Form("/group/nps/hhuang/analysis/DVCS_NPS2023/DVCS_analysis/pi0Calib_wf/calibTree/%s/prodTree_pass2_v3_%d_%d.root", Kine.Data(), runList[irun], iseg));
            // cout<<"Run:"<<runList[irun]<<"; Segment: "<<iseg<<endl;
        }

        //Set branches
        chain->SetBranchAddress("hms_Vx", &hms_Vx);
        chain->SetBranchAddress("hms_Vy", &hms_Vy);
        chain->SetBranchAddress("hms_Vz", &hms_Vz);
        chain->SetBranchAddress("caloev", &caloev);

        // Get kinematics variables________________________________________________
        Double_t Beam_energy = *db->GetEntry_d("BEAM_param_Energy", runList[irun]); // Beam energy in GeV
        Double_t HMS_mom = *db->GetEntry_d("SIMU_param_HMSmomentum", runList[irun]); // HMS central momentum in GeV/c
        Double_t HMS_angle = *db->GetEntry_d("SIMU_param_HMSangle", runList[irun]); // HMS angle in rad
        Double_t Target_amu = *db->GetEntry_d("TARGET_param_Amu", runList[irun]); // Target amu

        // The angle and distance of calorimeter
        Double_t NPS_dist = *db->GetEntry_d("CALO_geom_Dist", runList[irun]); // NPS distence in cm
        Double_t NPS_angle = *db->GetEntry_d("CALO_geom_Yaw", runList[irun]); // NPS angle in rad
        Int_t *caloMaskBlock = new Int_t[1080]; // Get the mask block information in NPS numbering Scheme
        caloMaskBlock = db->GetEntry_i("CALO_flag_MaskBlock", runList[irun]);

        // Elastic coefficients
        Double_t *coefElas = new Double_t[1080]; // Get the elastic coefficients (GeV/mV) in NPS numbering Scheme
        coefElas = db->GetEntry_d("CALO_calib_ElasCoef", runList[irun]);
        for(int iblk = 0; iblk < 1080; iblk++){
            if(coefElas[iblk] <= 0 && !caloMaskBlock[iblk]) coefElas[iblk] = 0.014; // if the coefficient <= 0, set it to 0.014 GeV/mV
        }

        nevt = chain->GetEntries();
        // Event loop: check the mass peak with old coefficients
        for(Int_t i = 0; i < nevt; i++){ // Event loop: check the mass peak with old coefficients
            chain->GetEntry(i);
            if (i % 1000 == 0) cout << i << "/" << nevt << endl;

            //Redo the cluster to update the pi0 mass 
            for(int iblk = 0; iblk < 1080; iblk++) hasBlock[iblk] = false;

            Int_t nCaloBlock = caloev->GetNbBlocks();
            for(int iblk = 0; iblk < nCaloBlock; iblk++){
                TCaloBlock *block = caloev->GetBlock(iblk);
                Int_t nb = block->GetBlockNumber(); // Block number in simulation numbering scheme
                Float_t energy = block->GetEnergy(0)*pi0_coef[bnConv_OldToNew(nb)]/coefElas[bnConv_OldToNew(nb)]; // Replace the coefficient

                block->Erase("");
                block->AddPulse(energy, 0);
                hasBlock[nb] = true;
            }
            
            for(Int_t iblk = 0; iblk < 1080; iblk++){
                if(hasBlock[iblk] == false) caloev->AddBlock(iblk);
            }

            // caloev->TriggerSim(0.5);
            caloev->DoClustering(-3,3); // Time window (-1.5,1.5) [ns]
            Int_t nclus = caloev->GetNbClusters();
            // cout<<nclus<<endl;

            for (Int_t iclus = 0; iclus < nclus; iclus++){
                caloev->GetCluster(iclus)->Analyze(1, -3, 3);
            }

            if(nclus != 2) continue;

            ev->SetVertex(hms_Vx, hms_Vy, hms_Vz);
            ev->GetGeometry()->SetCaloTheta(-1*NPS_angle);
            ev->GetGeometry()->SetCaloDist(NPS_dist);
            ev->SetCaloEvent(caloev);

            // Look for seed blocks for cut
            Double_t maxene;
            Int_t seed1, seed2;
            for(int iclus = 0; iclus < 2; iclus++){
                maxene = -1;
                if(iclus == 0) seed1 = -1;
                else seed2 = -1;
                for (Int_t k = 0; k < caloev->GetCluster(iclus)->GetClusSize(); k++){
                    TCaloBlock *clusblock = caloev->GetCluster(iclus)->GetBlock(k);
                    Double_t blockene = clusblock->GetBlockEnergy();
                    if (blockene > maxene){
                        maxene = blockene;
                        if(iclus == 0) seed1 = clusblock->GetBlockNumber();
                        else seed2 = clusblock->GetBlockNumber();
                    }
                }
            }

            // Check if the seeds are reasonable
            if(seed1 <  0 || seed2 < 0 || seed1 > 1079 || seed2 > 1079){
                cout<<"ERROR: Got a wrong seed block!!!"<<endl;
                return;
            }

            if(!(passAccCut(seed1) && passAccCut(seed2))) continue;

            TLorentzVector photon1 = ev->GetPhoton(0, 7, 0);
            TLorentzVector photon2 = ev->GetPhoton(1, 7, 0);
            Double_t m_pi0_recons = (photon1 + photon2).M();

            if(photon1.E() > 1.4 && photon2.E() > 1.4) h_pi0M->Fill(m_pi0_recons);

            ev->Reset();
        } // loop of events for old mass plots 

        f_fit->SetParLimits(1, 0.12, 0.15);
        f_fit->SetParLimits(2, 0.00001, 0.03);
        h_pi0M->Fit(f_fit, "R");

        N = f_fit->GetParameter(0)/h_pi0M->GetBinWidth(1);
        N_err = f_fit->GetParError(0)/h_pi0M->GetBinWidth(1);
        mean = f_fit->GetParameter(1);
        mean_err = f_fit->GetParError(1);
        sigm = f_fit->GetParameter(2);
        sigm_err = f_fit->GetParError(2);

        Int_t ibin = runList[irun]-first_run+1;
        h_pi0MeanvsRun_old->SetBinContent(ibin, mean*1000);
        h_pi0MeanvsRun_old->SetBinError(ibin, mean_err*1000);
        h_pi0SigmvsRun_old->SetBinContent(ibin, sigm*1000);
        h_pi0SigmvsRun_old->SetBinError(ibin, sigm_err*1000);

        outfile->cd();
        h_pi0M->Write(Form("h_pi0M_old_%d", runList[irun]));
        h_pi0M->Reset();

        // Update new coefficients for this run
        Double_t pi0_coef_new[1080];
        for(int iblk = 0; iblk < 1080; iblk++) pi0_coef_new[iblk] = 0;
        for(int iblk = 0; iblk < 1080; iblk++) pi0_coef_new[iblk] = pi0_coef[iblk]*(m_pi0 / mean); // in NPS numbering scheme

        // output to text file
        system(Form("rm -f Result/%s_pass2_v3/coef_pi0Calib_%d.txt", filename.Data(), runList[irun]));
        ofstream fcoef_out(Form("Result/%s_pass2_v3/coef_pi0Calib_%d.txt", filename.Data(), runList[irun]));
        for(int iblk = 0; iblk < 1080; iblk++) fcoef_out<<pi0_coef_new[iblk]<<endl;
        fcoef_out.close();

        // output to text file for nps replay
        system(Form("rm -f Result/%s_pass2_v3/coef_pi0Calib_%d_npsreplay.txt", filename.Data(), runList[irun]));
        ofstream fout_coef_hcana(Form("Result/%s_pass2_v3/coef_pi0Calib_%d_npsreplay.txt", filename.Data(), runList[irun]));
        fout_coef_hcana<<"nps_cal_arr_gain_cor = ";
        for(int iblk = 0; iblk < 1080; iblk++){
            if((iblk+1) % 30 == 0) fout_coef_hcana<<pi0_coef_new[iblk]<<", \n           ";
            else fout_coef_hcana<<pi0_coef_new[iblk]<<", ";
        }

        // Event loop: check the mass peak with updated coefficients
        for(Int_t i = 0; i < nevt; i++){ // Event loop: check the mass peak with updated coefficients
            chain->GetEntry(i);
            if (i % 1000 == 0) cout << i << "/" << nevt << endl;

            //Redo the cluster to update the pi0 mass 
            for(int iblk = 0; iblk < 1080; iblk++) hasBlock[iblk] = false;

            Int_t nCaloBlock = caloev->GetNbBlocks();
            for(int iblk = 0; iblk < nCaloBlock; iblk++){
                TCaloBlock *block = caloev->GetBlock(iblk);
                Int_t nb = block->GetBlockNumber(); // Block number in simulation numbering scheme
                Float_t energy = block->GetEnergy(0)*pi0_coef_new[bnConv_OldToNew(nb)]/coefElas[bnConv_OldToNew(nb)]; // Replace the coefficient

                block->Erase("");
                block->AddPulse(energy, 0);
                hasBlock[nb] = true;
            }
            
            for(Int_t iblk = 0; iblk < 1080; iblk++){
                if(hasBlock[iblk] == false) caloev->AddBlock(iblk);
            }

            // caloev->TriggerSim(0.5);
            caloev->DoClustering(-3,3); // Time window (-1.5,1.5) [ns]
            Int_t nclus = caloev->GetNbClusters();
            // cout<<nclus<<endl;

            for (Int_t iclus = 0; iclus < nclus; iclus++){
                caloev->GetCluster(iclus)->Analyze(1, -3, 3);
            }

            if(nclus != 2) continue;

            ev->SetVertex(hms_Vx, hms_Vy, hms_Vz);
            ev->GetGeometry()->SetCaloTheta(-1*NPS_angle);
            ev->GetGeometry()->SetCaloDist(NPS_dist);
            ev->SetCaloEvent(caloev);

            // Look for seed blocks for cut
            Double_t maxene;
            Int_t seed1, seed2;
            for(int iclus = 0; iclus < 2; iclus++){
                maxene = -1;
                if(iclus == 0) seed1 = -1;
                else seed2 = -1;
                for (Int_t k = 0; k < caloev->GetCluster(iclus)->GetClusSize(); k++){
                    TCaloBlock *clusblock = caloev->GetCluster(iclus)->GetBlock(k);
                    Double_t blockene = clusblock->GetBlockEnergy();
                    if (blockene > maxene){
                        maxene = blockene;
                        if(iclus == 0) seed1 = clusblock->GetBlockNumber();
                        else seed2 = clusblock->GetBlockNumber();
                    }
                }
            }

            // Check if the seeds are reasonable
            if(seed1 <  0 || seed2 < 0 || seed1 > 1079 || seed2 > 1079){
                cout<<"ERROR: Got a wrong seed block!!!"<<endl;
                return;
            }

            if(!(passAccCut(seed1) && passAccCut(seed2))) continue;

            TLorentzVector photon1 = ev->GetPhoton(0, 7, 0);
            TLorentzVector photon2 = ev->GetPhoton(1, 7, 0);
            Double_t m_pi0_recons = (photon1 + photon2).M();

            if(photon1.E() > 1.4 && photon2.E() > 1.4) h_pi0M->Fill(m_pi0_recons);

            ev->Reset();
        } // loop of events updated mass plot

        f_fit->SetParLimits(1, 0.12, 0.15);
        f_fit->SetParLimits(2, 0.00001, 0.03);
        h_pi0M->Fit(f_fit, "R");

        N = f_fit->GetParameter(0)/h_pi0M->GetBinWidth(1);
        N_err = f_fit->GetParError(0)/h_pi0M->GetBinWidth(1);
        mean = f_fit->GetParameter(1);
        mean_err = f_fit->GetParError(1);
        sigm = f_fit->GetParameter(2);
        sigm_err = f_fit->GetParError(2);

        h_pi0MeanvsRun_new->SetBinContent(ibin, mean*1000);
        h_pi0MeanvsRun_new->SetBinError(ibin, mean_err*1000);
        h_pi0SigmvsRun_new->SetBinContent(ibin, sigm*1000);
        h_pi0SigmvsRun_new->SetBinError(ibin, sigm_err*1000);

        outfile->cd();
        h_pi0M->Write(Form("h_pi0M_new_%d", runList[irun]));
        h_pi0M->Reset();

        chain->Reset();
    } // loop of runs

    outfile->cd();
    h_pi0MeanvsRun_old->Write();
    h_pi0MeanvsRun_new->Write();
    h_pi0SigmvsRun_old->Write();
    h_pi0SigmvsRun_new->Write();
    outfile->Close();
}

bool passAccCut(int sim_blk_number)
{   
    bool pass = false;
    // Column and row numbers for the cut on seed blocks
    // these col and row numbers are in sim. numbering scheme
    Int_t sim_irow = sim_blk_number%nrow;
    Int_t sim_icol = (sim_blk_number-sim_irow)/nrow;
    // cout<<"irow: "<<sim_irow<<"; icol: "<<sim_icol<<endl;

    if(0 < sim_irow && sim_irow < 35 && 0 < sim_icol && sim_icol < ncol_inUse-1) pass = true;

    return pass;
}