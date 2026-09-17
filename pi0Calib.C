// Execute with ./pi0Calib.sh <kinematics name> <target flag> <iteration to start>
// see /group/nps/hhuang/analysis/DVCS_NPS2023/DVCS_analysis/MakeProdList/<kinematics name>_<target flag>_ProdList.txt
// Target flag: 0 for LH2, 1 for LD2, -1 for both

#include "/group/nps/hhuang/analysis/MyHeader/Analysis.h"
#include "/group/nps/hhuang/analysis/MyHeader/MyDB.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include "TVectorD.h"
#include "TProcessID.h"
#include <unistd.h>  // for getpid()

Double_t photEcut = 1.5; // every 0.1 GeV within 0.5-1.5 GeV, photon energy cut for pi0 calibration
// Get kinematics variables from the database
Double_t Beam_energy, HMS_mom, HMS_angle, Target_amu; // GeV, GeV, rad
Double_t NPS_dist, NPS_angle; // unit in cm and rad in DB
Double_t *coefElas = new Double_t[nblk]; // Elastic coefficients in NPS numbering scheme
Double_t xnut, ynut; // center of the NPS
void SetKineVariables(TDVCSDB *db, int run_number);

//Get seed blocks
int GetSeedBlocks(Double_t xclus, Double_t yclus);

// Check if the block is in the range for calibration (with sim. numbering scheme)
bool PassAccCut(const vector<int> &maskBlkList, const double &NPS_dist, int sim_blk_number);

// pi0 calibration function
void pi0Calib(TString Kine = "60_4b", int TargetFlag = 0, int startcycle = 0, int group = 1, int iIter = 1)
{   
    TH1::SetDefaultSumw2();
    const Int_t tree_idx = 10*(photEcut-0.5);

    gSystem->Load("/group/nps/hhuang/software/NPS_SOFT/libDVCS.so");
    TDVCSDB *db = new TDVCSDB("dvcs", "clrlpc", 3306, "hhuang", "");

    // Get center positions of NPS for acceptance cut, should be run-independent
    xnut = *db->GetEntry_d("CALO_geom_X0", 1500);
    ynut = *db->GetEntry_d("CALO_geom_Y0", 1500);

    // Set target information
    TString Tar;
    Double_t M_targ = -1;
    if(TargetFlag == 0){ Tar = "LH2"; M_targ = 0.938272; }
    else if(TargetFlag == 1){ Tar = "LD2"; M_targ = 0.939565; }
    else{
        cout<<"ERROR: Unknown target!! Please specify the target for calibration!!"<<endl;
        return;
    }

    // Make folder for the new calibration results
    TString prefix_1 = "wf";
    // TString prefix_1 = "pass2";
    TString prefix_2 = Form("cycle%d", startcycle); for(int ic = startcycle+1; ic < startcycle+group; ic++) prefix_2+=Form("_cycle%d", ic);
    TString filename = Form("x%s_%s_%s_%s", Kine.Data(), Tar.Data(), prefix_1.Data(), prefix_2.Data());
    if(iIter == 1) system(Form("mkdir -p Result/%s", filename.Data()));
    
    // Get the list of runs for calibration
    vector<int> runList;
    vector<int> nSegList;
    for(int ic = startcycle; ic < startcycle+group; ic++){
        TString runlistDir = "/group/nps/hhuang/analysis/DVCS_NPS2023/DVCS_analysis/pi0Calib_wf/calibList/ListPerCycle";
        TString listname = Form("x%s_%d_cycle%d.txt", Kine.Data(), TargetFlag, ic);
        ifstream fRunList(Form("%s/%s", runlistDir.Data(), listname.Data()));
        
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
        fRunList.close();
    }
    
    Int_t nRun = runList.size();
    Int_t first_run = runList[0];
    Int_t last_run = runList[nRun-1];

    // for(int irun = 0; irun < nRun; irun++){
    //     cout<<"Run: "<<runList[irun]<<"; nSeg: "<<nSegList[irun]<<endl;
    // }

    // Add the rootfiles to TChain
    TString calibTreeDir = "/group/nps/hhuang/analysis/DVCS_NPS2023/DVCS_analysis/pi0Calib_wf/calibTree";
    TChain *chain = new TChain(Form("t_prod_%d", tree_idx));
    for(int irun = 0; irun < nRun; irun++){
        for(int iseg = 0; iseg < nSegList[irun]; iseg++){
            // cout<<"Run:"<<runList[irun]<<"; Segment: "<<iseg<<endl;
            chain->Add(Form("%s/x%s/prodTree_%s_%d_%d.root", calibTreeDir.Data(), Kine.Data(), prefix_1.Data(), runList[irun], iseg));
            // cout<<"Run:"<<runList[irun]<<"; Segment: "<<iseg<<endl;
        }
    }

    // Get branches in the rootfiles
    Double_t runNb;
    Double_t H_react_x;
    Double_t H_react_y;
    Double_t H_react_z;
    Double_t H_gtr_px;
    Double_t H_gtr_py;
    Double_t H_gtr_pz;
    // Double_t m;
    // Double_t mm2;

    TCaloEvent *caloev = new TCaloEvent(first_run);
    chain->SetBranchAddress("g.runnum", &runNb);
    chain->SetBranchAddress("H_react_x", &H_react_x);
    chain->SetBranchAddress("H_react_y", &H_react_y);
    chain->SetBranchAddress("H_react_z", &H_react_z);
    // chain->SetBranchAddress("m", &m);
    // chain->SetBranchAddress("mm2", &mm2);
    chain->SetBranchAddress("caloev", &caloev);

    // Get kinematics variables
    SetKineVariables(db, first_run);

    // List of block to mask out (all blocks found in the run list)
    vector <int> maskBlkList; // block numbers to mask out in NPS numbering scheme
    Int_t *caloMaskBlock = new Int_t[nblk]; // Get the mask block information in NPS numbering Scheme
    for(int irun = 0; irun < nRun; irun++){
        caloMaskBlock = db->GetEntry_i("CALO_flag_MaskBlock", runList[irun]);
        for(int iblk = 0; iblk < nblk; iblk++){
            if(caloMaskBlock[iblk]) maskBlkList.push_back(iblk);
        }
    }
    sort(maskBlkList.begin(), maskBlkList.end()); // sort the block numbers
    auto last = unique(maskBlkList.begin(), maskBlkList.end()); // move the duplicate numbers to the end
    maskBlkList.erase(last, maskBlkList.end());  // remove the duplicate numbers

    // Print out some basic information for recording
    system(Form("rm -f Result/%s/calib_Info.txt", filename.Data())); // avoid appending old file
    ofstream info_Stream(Form("Result/%s/calib_Info.txt", filename.Data()));

    info_Stream<<"This is the calibration for KinC_x"<<Kine.Data()<<endl;
    info_Stream<<"Beam energy: "<<Beam_energy<<endl;
    info_Stream<<"Target in use: "<<GetTarName(Target_amu)<<endl;
    info_Stream<<" "<<endl;
    info_Stream<<"HMS momentum: "<<HMS_mom<<" GeV/c"<<endl;
    info_Stream<<"HMS angle: "<<HMS_angle*TMath::RadToDeg()<<" deg."<<endl;
    info_Stream<<" "<<endl;
    info_Stream<<"SHMS angle: "<<NPS_angle*TMath::RadToDeg()+16.3<<" deg."<<endl;
    info_Stream<<"NPS angle: "<<NPS_angle*TMath::RadToDeg()<<" deg."<<endl;
    info_Stream<<"NPS distance: "<<NPS_dist/100<<" m"<<endl;
    info_Stream<<" "<<endl;
    info_Stream<<"Photon energy cut: "<<photEcut<<" GeV"<<endl;
    info_Stream<<" "<<endl;
    info_Stream<<"Run list: calibList/"<<Kine.Data()<<"/"<<Tar.Data()<<".txt"<<endl;
    for(int irun = 0; irun < nRun; irun++) info_Stream<<"Run "<<runList[irun]<<endl;
    info_Stream<<" "<<endl;

    // Initial settings for calibration
    TDVCSEvent *ev = new TDVCSEvent(first_run);
    Bool_t hasBlock[nblk];

    // matrix and vectors used to compute epsilon
    Float_t B = 0;
    Float_t B2 = 0; // to calculate F for debugging
    Float_t LCD = 0;
    Float_t LCL = 0;
    Float_t lambda = 0;
    Float_t m_pi0_recons = 0;

    TMatrixD C_mat(nblk, nblk);
    TMatrixD C_mat_inv(nblk, nblk);
    TMatrixD eigen_vectors(nblk, nblk);
    TVectorD eigen_vals(nblk);
    Double_t eigen_vectors_norm[nblk];

    Float_t corr_pi0[nblk]; // correction factor from pi0 calibration in simulation numbering scheme
    Float_t block_energy_ratio[nblk];
    Float_t dmdepsilon[nblk];
    Float_t epsilon[nblk];
    for(Int_t i = 0; i < nblk; i++){
        block_energy_ratio[i] = 0;
        dmdepsilon[i] = 0;
        epsilon[i] = 0;
    }

    Float_t D[nblk];
    Float_t L[nblk];
    for(Int_t i = 0; i < nblk; i++){
        D[i] = 0;
        L[i] = 0;
        eigen_vals(i) = 0.0;
        eigen_vectors_norm[i] = 0.0;

        for(Int_t j = 0; j < nblk; j++){
            C_mat(i, j) = 0;
            C_mat_inv(i, j) = 0;
            eigen_vectors(i, j) = 0.0;
        }
    }

    // pi0 selection and fit its mass distribution (will be updated later)
    Float_t pi0Mass_min = 0.08;
    Float_t pi0Mass_max = 0.15;

    // Ouput files of each iterations for recording
    system(Form("rm -f Result/%s/log_Iteration_%d.txt", filename.Data(), iIter)); // avoid appending old file
    ofstream log_Stream(Form("Result/%s/log_Iteration_%d.txt", filename.Data(), iIter));
    
    log_Stream<<"++++++++++ This is iteration "<<iIter<<" of calibration ++++++++++"<<endl;

    Int_t nevt = chain->GetEntries();

    if(iIter == 1){ // initialize the correction factor and generate the mass distribution from clusters
        for(Int_t i = 0; i < nblk; i++) corr_pi0[i] = 1; // initialize the correction factor

        // Output histogram of invariant mass distribution
        TFile *output_iter0 = new TFile(Form("Result/%s/f_Mgg_%d.root", filename.Data(), iIter-1), "recreate");
        TH1F *h_pi0M_iter0 = new TH1F(Form("h_pi0M_%d", iIter-1), Form("M_{#gamma#gamma} (E_{#gamma} > %.1f GeV);M_{#gamma#gamma} [GeV];Counts", photEcut), 150, 0.05, 0.2);
        h_pi0M_iter0->Sumw2();

        Double_t last_runNb = first_run; // track the change of run number
        for(Int_t i = 0; i < nevt; i++){ // Event loop: get the mass distribution from clusters
            Int_t ObjectNumber = TProcessID::GetObjectCount();
            chain->GetEntry(i);
            if (i % 10000 == 0) cout << i << "/" << nevt << endl;

            // Update calorimeter geometry if the run number changes
            if(runNb != last_runNb){
                info_Stream<<"Run number changed from "<<last_runNb<<" to "<<runNb<<endl;
                info_Stream<<"Update calo distance and angle"<<endl;
                last_runNb = runNb; // update the last run number

                // The angle and distance of calorimeter
                NPS_dist = *db->GetEntry_d("CALO_geom_Dist", last_runNb); // NPS distence in cm
                NPS_angle = *db->GetEntry_d("CALO_geom_Yaw", last_runNb); // NPS angle in rad

                info_Stream<<"NPS angle: "<<NPS_angle*TMath::RadToDeg()<<" deg."<<endl;
                info_Stream<<"NPS distance: "<<NPS_dist/100<<" m"<<endl;
            }

            if(!(TMath::Abs(H_react_z) < 15)) continue;
            ev->SetVertex(H_react_x, H_react_y, H_react_z);
            ev->GetGeometry()->SetCaloTheta(-1*NPS_angle);
            ev->GetGeometry()->SetCaloDist(NPS_dist);
            ev->SetCaloEvent(caloev);

            if(caloev->GetNbClusters() != 2) continue; // Only look at two clusters events

            // Check seed blocks and apply acceptance cut
            Int_t seed1 = GetSeedBlocks(caloev->GetCluster(0)->GetX(), caloev->GetCluster(0)->GetY());
            Int_t seed2 = GetSeedBlocks(caloev->GetCluster(1)->GetX(), caloev->GetCluster(1)->GetY());
            if(!(PassAccCut(maskBlkList, NPS_dist, seed1) && PassAccCut(maskBlkList, NPS_dist, seed2))) continue;

            // Get reconstructed photons
            TLorentzVector photon1 = ev->GetPhoton(0, 7, 0);
            TLorentzVector photon2 = ev->GetPhoton(1, 7, 0);
            m_pi0_recons = (photon1 + photon2).M();
            if(photon1.E() > photEcut && photon2.E() > photEcut) h_pi0M_iter0->Fill(m_pi0_recons);

            ev->Reset();
            TProcessID::SetObjectCount(ObjectNumber);
        }// event loop

        Double_t meanfit_iter0 = h_pi0M_iter0->GetBinCenter(h_pi0M_iter0->GetMaximumBin());
        Double_t minfit_iter0 = meanfit_iter0 - 0.035; // or 0.035 if using a background function
        Double_t maxfit_iter0 = meanfit_iter0 + 0.035;
        TF1 *f_fit_iter0 = new TF1(Form("f_fit_%d", iIter-1), "gausn(0)+pol1(3)", minfit_iter0, maxfit_iter0);
        // f_fit_iter0->SetParLimits(0, 100, 100000);
        f_fit_iter0->SetParLimits(1, 0.1, 0.15);
        f_fit_iter0->SetParLimits(2, 0.00001, 0.03);
        h_pi0M_iter0->Fit(f_fit_iter0, "R");

        Float_t N = f_fit_iter0->GetParameter(0)/h_pi0M_iter0->GetBinWidth(1);
        Float_t N_err = f_fit_iter0->GetParError(0)/h_pi0M_iter0->GetBinWidth(1);
        Float_t mean = f_fit_iter0->GetParameter(1);
        Float_t mean_err = f_fit_iter0->GetParError(1);
        Float_t sigm = f_fit_iter0->GetParameter(2);
        Float_t sigm_err = f_fit_iter0->GetParError(2);
        
        // update the pi0 mass window for the calibration
        pi0Mass_min = mean-3*sigm;
        pi0Mass_max = mean+3*sigm;

        // Save histogram and fitting function
        output_iter0->cd();
        h_pi0M_iter0->Write();
        f_fit_iter0->Write();
        output_iter0->Close();

        log_Stream<<"The pi0 mass window before calibration is: "<<pi0Mass_min<<" to "<<pi0Mass_max<<" GeV"<<endl;
        log_Stream<<"The estmated number of pi0 = "<<N<<" \u00B1 "<<N_err<<endl;
        log_Stream<<"Mean = "<<mean<<" \u00B1 "<<mean_err<<" GeV"<<endl;
        log_Stream<<"Sigma = "<<sigm<<" \u00B1 "<<sigm_err<<" GeV"<<endl;
    }
     
    else{// Get the correction factor and mass distrubution from the previous iteration when iIter > 1
        ifstream fcorr_old(Form("Result/%s/corr_pi0.txt", filename.Data()));
        for(int iblk = 0; iblk < nblk; iblk++) fcorr_old>>corr_pi0[iblk];

        // Get the mass distribution with the histogram from previous iteration
        TFile *input = TFile::Open(Form("Result/%s/f_Mgg_%d.root", filename.Data(), iIter-1));
        TH1F *h_pi0M_last = (TH1F*)input->Get(Form("h_pi0M_%d", iIter-1));

        Double_t meanfit_last = h_pi0M_last->GetBinCenter(h_pi0M_last->GetMaximumBin());
        Double_t minfit_last = meanfit_last - 0.035;
        Double_t maxfit_last = meanfit_last + 0.035;
        TF1 *f_fit_last = new TF1(Form("f_fit_%d", iIter-1), "gausn(0)+pol1(3)", minfit_last, maxfit_last);
        // f_fit_last->SetParLimits(0, 100, 100000);
        f_fit_last->SetParLimits(1, 0.1, 0.15);
        f_fit_last->SetParLimits(2, 0.00001, 0.03);
        h_pi0M_last->Fit(f_fit_last, "RN");

        Float_t N = f_fit_last->GetParameter(0)/h_pi0M_last->GetBinWidth(1);
        Float_t N_err = f_fit_last->GetParError(0)/h_pi0M_last->GetBinWidth(1);
        Float_t mean = f_fit_last->GetParameter(1);
        Float_t mean_err = f_fit_last->GetParError(1);
        Float_t sigm = f_fit_last->GetParameter(2);
        Float_t sigm_err = f_fit_last->GetParError(2);
        
        // update the pi0 mass window for the calibration
        pi0Mass_min = mean-3*sigm;
        pi0Mass_max = mean+3*sigm;

        log_Stream<<"The pi0 mass window after last iteration is: "<<pi0Mass_min<<" to "<<pi0Mass_max<<" GeV"<<endl;
        log_Stream<<"The estmated number of pi0 = "<<N<<" \u00B1 "<<N_err<<endl;
        log_Stream<<"Mean = "<<mean<<" \u00B1 "<<mean_err<<" GeV"<<endl;
        log_Stream<<"Sigma = "<<sigm<<" \u00B1 "<<sigm_err<<" GeV"<<endl;
        log_Stream<<"+++++++++++++++++++++++++"<<endl;
        log_Stream<<" "<<endl;
    }

    if(iIter <= 8){
        log_Stream<<"Start the pi0 calibration......"<<endl;
        Double_t last_runNb = first_run; // track the change of run number
        for (Int_t i = 0; i < nevt; i++){ // Event loop: pi0 calibration
            Int_t ObjectNumber = TProcessID::GetObjectCount();
            chain->GetEntry(i);
            if (i % 10000 == 0) cout << i << "/" << nevt << endl;

            // Update calorimeter geometry if the run number changes
            if(runNb != last_runNb){
                info_Stream<<"Run number changed from "<<last_runNb<<" to "<<runNb<<endl;
                info_Stream<<"Update calo distance and angle"<<endl;
                last_runNb = runNb; // update the last run number

                // The angle and distance of calorimeter
                NPS_dist = *db->GetEntry_d("CALO_geom_Dist", last_runNb); // NPS distence in cm
                NPS_angle = *db->GetEntry_d("CALO_geom_Yaw", last_runNb); // NPS angle in rad

                info_Stream<<"NPS angle: "<<NPS_angle*TMath::RadToDeg()<<" deg."<<endl;
                info_Stream<<"NPS distance: "<<NPS_dist/100<<" m"<<endl;
            }

            //Redo the clustering to update the cluster energy and pi0 mass with correction factors
            for(int iblk = 0; iblk < nblk; iblk++) hasBlock[iblk] = false;

            Int_t nCaloBlock = caloev->GetNbBlocks();
            for(int iblk = 0; iblk < nCaloBlock; iblk++){
                TCaloBlock *block = caloev->GetBlock(iblk);
                Int_t nb = block->GetBlockNumber();
                Float_t energy = block->GetEnergy(0)*corr_pi0[nb];

                block->Erase("");
                block->AddPulse(energy, 0);
                hasBlock[nb] = true;
            }
            
            for(Int_t iblk = 0; iblk < nblk; iblk++){
                if(hasBlock[iblk] == false) caloev->AddBlock(iblk);
            }

            // caloev->TriggerSim(0.2);
            caloev->DoClustering(-3,3); // Time window (-1.5,1.5) [ns]
            Int_t nclus = caloev->GetNbClusters();
            // cout<<nclus<<endl;

            for (Int_t iclus = 0; iclus < nclus; iclus++){
                caloev->GetCluster(iclus)->Analyze(1, -3, 3);
            }

            // 2 clusters and photon energy cut
            if(nclus != 2) continue;
            if(!(caloev->GetCluster(0)->GetEnergy() > photEcut && caloev->GetCluster(1)->GetEnergy() > photEcut)) continue;

            // Check seed blocks and apply acceptance cut
            Int_t seed1 = GetSeedBlocks(caloev->GetCluster(0)->GetX(), caloev->GetCluster(0)->GetY());
            Int_t seed2 = GetSeedBlocks(caloev->GetCluster(1)->GetX(), caloev->GetCluster(1)->GetY());
            if(!(PassAccCut(maskBlkList, NPS_dist, seed1) && PassAccCut(maskBlkList, NPS_dist, seed2))) continue;

            if(!(TMath::Abs(H_react_z) < 15)) continue;
            ev->SetVertex(H_react_x, H_react_y, H_react_z);
            ev->GetGeometry()->SetCaloTheta(-1*NPS_angle);
            ev->GetGeometry()->SetCaloDist(NPS_dist);
            ev->SetCaloEvent(caloev);

            // Get reconstructed photons
            TLorentzVector photon1 = ev->GetPhoton(0, 7, 0);
            TLorentzVector photon2 = ev->GetPhoton(1, 7, 0);
            m_pi0_recons = (photon1 + photon2).M();
            if(!(pi0Mass_min < m_pi0_recons && m_pi0_recons < pi0Mass_max)) continue;

            // caloev->GetCluster(0)->Print();
            for (Int_t k = 0; k < caloev->GetCluster(0)->GetClusSize(); k++){
                TCaloBlock *clusblock = caloev->GetCluster(0)->GetBlock(k);
                Float_t blockene = clusblock->GetBlockEnergy();
                // cout<<clusblock->GetBlockNumber()<<" == "<<clusblock->GetBlockEnergy()<<endl;

                // get block_energy_ratio from cluster 1
                block_energy_ratio[clusblock->GetBlockNumber()] = blockene/caloev->GetCluster(0)->GetEnergy();
            }

            for (Int_t k = 0; k < caloev->GetCluster(1)->GetClusSize(); k++){
                TCaloBlock *clusblock = caloev->GetCluster(1)->GetBlock(k);
                Float_t blockene = clusblock->GetBlockEnergy();

                // get block_energy_ratio from cluster 2
                block_energy_ratio[clusblock->GetBlockNumber()] = blockene/caloev->GetCluster(1)->GetEnergy();
            }

            // Check if the ration to cluster energy are reasonable
            bool hasBug = false;
            for(int k = 0; k < nblk; k++){
                if(block_energy_ratio[k] < 0 || block_energy_ratio[k] > 1){
                    hasBug = true;
                    cout<<"ERROR: There is a bug with block energy ratio of block "<<k<<", Skip this event."<<endl;
                }
            }
            if(hasBug) continue;
            
            // Calibration: compute B
            B += m_pi0_recons*m_pi0_recons - m_pi0*m_pi0;
            B2 += TMath::Power(m_pi0_recons*m_pi0_recons - m_pi0*m_pi0, 2);

            // compute C_mat matrix and vectors D and L
            for (Int_t j = 0; j < nblk; j++){
                dmdepsilon[j] = m_pi0_recons * m_pi0_recons * block_energy_ratio[j];
            }

            for(Int_t j = 0; j < nblk; j++){
                if (dmdepsilon[j] != 0){
                    D[j] -= (m_pi0_recons * m_pi0_recons - m_pi0 * m_pi0) * dmdepsilon[j];
                    L[j] += dmdepsilon[j];

                    for (Int_t k = 0; k < nblk; k++){
                        C_mat(j, k) += dmdepsilon[j] * dmdepsilon[k];
                    } // end for k
                }// end if
            }// end for j

            // reinitialization for next event
            ev->Reset();
            m_pi0_recons = 0.0;
            
            for(Int_t j = 0; j < nblk; j++){
                block_energy_ratio[j] = 0.0;
                dmdepsilon[j] = 0.0;
            }
            TProcessID::SetObjectCount(ObjectNumber);
        } // Event loop

        // Calculate calibration coefficients
        Int_t count_zero = 0;
        for(Int_t j = 0; j < nblk; j++){
            for(Int_t k = 0; k < nblk; k++){
                if(j != k) continue;
                if(C_mat(j, k) == 0){
                    cout<<"Warning: block "<<bnConv_OldToNew(j)<<" is 0 on the diagonal"<<endl;
                    count_zero += 1;
                }
            }
        }

        TMatrixD C_mat_clone(C_mat); // clone the origin matrix for checking the invertion

        if(count_zero > 0){
            //use spectral decomposition to invert matrix, and remove eigenvalues compatible with 0
            eigen_vectors=C_mat.EigenVectors(eigen_vals);

            //print eigen values of C_mat for debug and stuff
            // cout << "C_mat eigen values are:" << endl;
            // for(Int_t i=0;i<nblk;i++){
            //     cout << eigen_vals(i) << endl;
            // }

            //save the eigen values of C_mat in a txt file to keep as logs (for debug and stuff)
            // eigen_vals_stream << "iteration num: " << iter << "/" << iterations << " ; C_mat eigen values are:" << endl;
            // for(Int_t i=0;i<nblk;i++){
            //     eigen_vals_stream << eigen_vals(i) << endl;
            // }
            // eigen_vals_stream << endl;


            for(Int_t i=0;i<nblk;i++){
                for(Int_t j=0;j<nblk;j++){
                    eigen_vectors_norm[i]+=eigen_vectors(j,i)*eigen_vectors(j,i);//norm2, squared
                }
            }

            for(Int_t i=0;i<nblk;i++){
                for(Int_t j=0;j<nblk;j++){
                    for(Int_t k=0;k<nblk;k++){
                        if(eigen_vals(k)>0.01 || eigen_vals(k)<-0.01){//choice on the limit value of eigen values to be checked...
                            C_mat_inv(i,j)+=eigen_vectors(i,k)*eigen_vectors(j,k)/(eigen_vals(k)*sqrt(eigen_vectors_norm[k])*sqrt(eigen_vectors_norm[k]));
                        }
                    }
                }
            }
        }

        else C_mat_inv = C_mat.Invert(); // inversion

        TMatrixD C_check =  C_mat_inv*C_mat_clone; //debug
        // C_check.Print(); //debug
        
        // computation of epsilons
        for (Int_t i = 0; i < nblk; i++){
            for (Int_t j = 0; j < nblk; j++){
                LCD += L[i] * C_mat_inv(i, j) * D[j];
                LCL += L[i] * C_mat_inv(i, j) * L[j];
            }
        }

        lambda = (B + LCD) / (LCL);

        for (Int_t i = 0; i < nblk; i++){
            for (Int_t j = 0; j < nblk; j++){
                epsilon[i] += C_mat_inv(i, j) * (D[j] - lambda * L[j]);
            }
        }

        // update correction factors with epsilon
        for (Int_t i = 0; i < nblk; i++) corr_pi0[i] *= (1 + epsilon[i]);

        // update correction coefficients for those outside acceptance with the average of the others
        Double_t sum_inAcc = 0;
        Double_t count_inAcc = 0;
        for (Int_t i = 0; i < nblk; i++){
            if(PassAccCut(maskBlkList, NPS_dist, i)){
                sum_inAcc+=coefElas[bnConv_OldToNew(i)]*corr_pi0[i];
                count_inAcc+=1.;
            }
            // cout<<corr_pi0[i]<<" : "<<sum_inAcc<<endl;
        }

        for (Int_t i = 0; i < nblk; i++){
            if(!PassAccCut(maskBlkList, NPS_dist, i)) corr_pi0[i] = (sum_inAcc/count_inAcc)/coefElas[bnConv_OldToNew(i)];
        }

        // save the correction factors corr_pi0 
        system(Form("rm -f Result/%s/corr_pi0.txt", filename.Data())); // avoid appending old file
        ofstream corr_pi0_stream(Form("Result/%s/corr_pi0.txt", filename.Data()));
        for (Int_t i = 0; i < nblk; i++){
            corr_pi0_stream<<corr_pi0[i]<<endl;
            // cout << corr_pi0[i] << endl;
        }
        corr_pi0_stream.close();

        // Output for debugging____________________________________________________
        log_Stream<<endl;
        log_Stream<<"+++++ information for debugging +++++"<<endl;
        log_Stream<<"The physical mass of pi0 is set as "<<m_pi0<<" GeV"<<endl;
        log_Stream<<"LCD = "<<LCD<<endl;
        log_Stream<<"LCL = "<<LCL<<endl;
        log_Stream<<"B = "<<B<<endl;
        log_Stream<<"B2 = "<<B2<<endl;
        log_Stream<<"lambda = "<<lambda<<endl;
        log_Stream<<"F = "<<B2+2*lambda*B<<endl;
        log_Stream<<"+++++++++++++++++++++++++++++++++++++"<<endl;
        log_Stream<<endl;
        
        TH2D *hh_mat = new TH2D(C_mat_clone); // matrix
        TH2D *hh_mat_inv = new TH2D(C_mat_inv); // invert matrix
        TH2D *hh_mat_mult = new TH2D(C_check); // their multiplication
        TH1F *h_epsilon = new TH1F(Form("h_epsilon_%d", iIter),"#epsilon of #pi^{0} calibration;Block number in simulation numbering scheme;#epsilon", 1080, -0.5, 1079.5); // epsilon
        for (Int_t iblk = 0; iblk < nblk; iblk++) h_epsilon->SetBinContent(iblk+1, epsilon[iblk]);

        TFile *output_debug = new TFile(Form("Result/%s/f_Debug_%d.root", filename.Data(), iIter), "recreate");
        hh_mat->Write(Form("hh_mat_%d", iIter));
        hh_mat_inv->Write(Form("hh_mat_inv_%d", iIter));
        hh_mat_mult->Write(Form("hh_mat_mult_%d", iIter));
        h_epsilon->Write();
        output_debug->Close();
        // end of output for debugging_______________________________________________

        log_Stream<<"The pi0 calibration of iteration "<<iIter<<" finished"<<endl;
        log_Stream<<"Update the invariant mass distribution with last iteration......"<<endl;

        
        // Event loop: check the new mass position of pi0__________________________________________________________
        TFile *output = new TFile(Form("Result/%s/f_Mgg_%d.root", filename.Data(), iIter), "recreate");
        TH1F *h_pi0M = new TH1F(Form("h_pi0M_%d", iIter), Form("M_{#gamma#gamma} (E_{#gamma} > %.1f GeV);M_{#gamma#gamma} [GeV];Counts", photEcut), 150, 0.05, 0.2);
        h_pi0M->Sumw2();

        last_runNb = first_run; // track the change of run number
        for (Int_t i = 0; i < nevt; i++){ // Event loop: update the invariant mass distribution after the last iteration
            Int_t ObjectNumber = TProcessID::GetObjectCount();
            chain->GetEntry(i);
            if (i % 10000 == 0) cout << i << "/" << nevt << endl;

            // Update calorimeter geometry if the run number changes
            if(runNb != last_runNb){
                info_Stream<<"Run number changed from "<<last_runNb<<" to "<<runNb<<endl;
                info_Stream<<"Update calo distance and angle"<<endl;
                last_runNb = runNb; // update the last run number

                // The angle and distance of calorimeter
                NPS_dist = *db->GetEntry_d("CALO_geom_Dist", last_runNb); // NPS distence in cm
                NPS_angle = *db->GetEntry_d("CALO_geom_Yaw", last_runNb); // NPS angle in rad

                info_Stream<<"NPS angle: "<<NPS_angle*TMath::RadToDeg()<<" deg."<<endl;
                info_Stream<<"NPS distance: "<<NPS_dist/100<<" m"<<endl;
            }

            //Redo the cluster to update the pi0 mass 
            for(int iblk = 0; iblk < nblk; iblk++) hasBlock[iblk] = false;

            Int_t nCaloBlock = caloev->GetNbBlocks();
            for(int iblk = 0; iblk < nCaloBlock; iblk++){ 
                TCaloBlock *block = caloev->GetBlock(iblk);
                Int_t nb = block->GetBlockNumber(); // Here the block number is in simulation numbering scheme
                Float_t energy = block->GetEnergy(0)*corr_pi0[nb];

                block->Erase("");
                block->AddPulse(energy, 0);
                hasBlock[nb] = true;
            }
            
            for(Int_t iblk = 0; iblk < nblk; iblk++){
                if(hasBlock[iblk] == false) caloev->AddBlock(iblk);
            }

            // caloev->TriggerSim(0.2);
            caloev->DoClustering(-3,3); // Time window (-1.5,1.5) [ns]
            Int_t nclus = caloev->GetNbClusters();
            // cout<<nclus<<endl;

            for (Int_t iclus = 0; iclus < nclus; iclus++){
                caloev->GetCluster(iclus)->Analyze(1, -3, 3);
            }

            if(nclus != 2) continue;

            if(!(TMath::Abs(H_react_z) < 15)) continue;
            ev->SetVertex(H_react_x, H_react_y, H_react_z);
            ev->GetGeometry()->SetCaloTheta(-1*NPS_angle);
            ev->GetGeometry()->SetCaloDist(NPS_dist);
            ev->SetCaloEvent(caloev);

            // Check seed blocks and apply acceptance cut
            Int_t seed1 = GetSeedBlocks(caloev->GetCluster(0)->GetX(), caloev->GetCluster(0)->GetY());
            Int_t seed2 = GetSeedBlocks(caloev->GetCluster(1)->GetX(), caloev->GetCluster(1)->GetY());
            if(!(PassAccCut(maskBlkList, NPS_dist, seed1) && PassAccCut(maskBlkList, NPS_dist, seed2))) continue;

            // Get reconstructed photons
            TLorentzVector photon1 = ev->GetPhoton(0, 7, 0);
            TLorentzVector photon2 = ev->GetPhoton(1, 7, 0);
            m_pi0_recons = (photon1 + photon2).M();

            if(photon1.E() > photEcut && photon2.E() > photEcut) h_pi0M->Fill(m_pi0_recons);

            ev->Reset();
            TProcessID::SetObjectCount(ObjectNumber);
        }// event loop

        Double_t meanfit = h_pi0M->GetBinCenter(h_pi0M->GetMaximumBin());
        Double_t minfit = meanfit - 0.035;
        Double_t maxfit = meanfit + 0.035;
        TF1 *f_fit = new TF1(Form("f_fit_%d", iIter-1), "gausn(0)+pol1(3)", minfit, maxfit);
        // f_fit->SetParLimits(0, 100, 100000);
        f_fit->SetParLimits(1, 0.1, 0.15);
        f_fit->SetParLimits(2, 0.00001, 0.03);
        h_pi0M->Fit(f_fit, "R");

        Float_t N = f_fit->GetParameter(0)/h_pi0M->GetBinWidth(1);
        Float_t N_err = f_fit->GetParError(0)/h_pi0M->GetBinWidth(1);
        Float_t mean = f_fit->GetParameter(1);
        Float_t mean_err = f_fit->GetParError(1);
        Float_t sigm = f_fit->GetParameter(2);
        Float_t sigm_err = f_fit->GetParError(2);
        
        // update the pi0 mass window for the calibration
        pi0Mass_min = mean-3*sigm;
        pi0Mass_max = mean+3*sigm;

        // Save histogram and fitting function
        output->cd();
        h_pi0M->Write();
        output->Close();

        log_Stream<<"The pi0 mass window after calibration is: "<<pi0Mass_min<<" to "<<pi0Mass_max<<" GeV"<<endl;
        log_Stream<<"The estmated number of pi0 = "<<N<<" \u00B1 "<<N_err<<endl;
        log_Stream<<"Mean = "<<mean<<" \u00B1 "<<mean_err<<" GeV"<<endl;
        log_Stream<<"Sigma = "<<sigm<<" \u00B1 "<<sigm_err<<" GeV"<<endl;

        log_Stream<<"The new correction factors are:"<<endl;
        for (Int_t i = 0; i < nblk; i++){
            log_Stream<<corr_pi0[i]<<endl;;
        }
        // log_Stream<<"pi0 calibration completed"<<endl;
        log_Stream<<"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
    }

    if(iIter == 9){ // output the coefficients and make histogram after finishing the calibration
        TFile *outfile = new TFile(Form("Result/%s/CalibSummary.root", filename.Data()), "recreate");

        // Input correction factor
        ifstream fcorr_old(Form("Result/%s/corr_pi0.txt", filename.Data()));
        for(int iblk = 0; iblk < nblk; iblk++) fcorr_old>>corr_pi0[iblk]; // simulation numbering scheme
        // Convert the correction factor to NPS numbering scheme
        Double_t corr_pi0_nps[nblk];
        for(int iblk = 0; iblk < nblk; iblk++) corr_pi0_nps[iblk] = corr_pi0[bnConv_NewToOld(iblk)]; // convert to NPS numbering scheme

        // Calculate the actural coefficients after this calibration then output to text file
        Double_t pi0_coef[nblk];
        system(Form("rm -f Result/%s/coef_pi0Calib_temp.txt", filename.Data()));
        ofstream fout_coef(Form("Result/%s/coef_pi0Calib_temp.txt", filename.Data()));
        for(int iblk = 0; iblk < nblk; iblk++){ // in NPS numbering scheme
            pi0_coef[iblk] = coefElas[iblk]*corr_pi0_nps[iblk];
            fout_coef<<pi0_coef[iblk]<<endl;
        }
        
        // Make histograms for these coefficients
        TH1F *h_coEff = new TH1F("h_coEff", "Calibration coefficients;PMT number;coefficients", 1080, -0.5,  1079.5);
        TH2F *hh_coEff = new TH2F("hh_coEff", "Calibration coefficients; Column number; Row number; ", 30, -0.5, 29.5, 36, -0.5, 35.5);
        for(int iblk = 0; iblk < nblk; iblk++){
            Int_t icol = iblk%30;
            Int_t irow = (iblk-icol)/30;

            // Fill the calibration coefficients into histogram
            h_coEff->SetBinContent(iblk+1, pi0_coef[iblk]);
            hh_coEff->SetBinContent(icol+1, irow+1, pi0_coef[iblk]);
        }

        // Then output the coefficients for each run with a further correction

        nevt = -1; // Initialize
        chain->Reset();
        Beam_energy = -999;
        HMS_mom = -999;
        HMS_angle = -999;
        Target_amu = -999;
        NPS_dist = -999;
        NPS_angle = -999;
        maskBlkList.clear();
        for(int iblk = 0; iblk < nblk; iblk++){
            caloMaskBlock[iblk] = 1;
            coefElas[iblk] = 0;
        }
        Float_t N, N_err, mean, mean_err, sigm, sigm_err; // For fitting pi0 mass peak
        
        TH1F *h_pi0M_final = new TH1F("h_pi0M_final", Form("M_{#gamma#gamma} (E_{#gamma1} > %.1f GeV);M_{#gamma#gamma} [GeV];Counts", photEcut), 150, 0.05, 0.2);
        h_pi0M_final->Sumw2();

        TH1F *h_pi0Mx2_final = new TH1F("h_pi0Mx2_final", Form("M_{x}^{2} (E_{#gamma1} > %.1f GeV);M_{x}^{2} [GeV^{2}];Counts", photEcut), 100, 0, 1.5);
        h_pi0Mx2_final->Sumw2();

        TF1 *f_fit_final = new TF1("f_fit_final", "gausn(0)+pol1(3)", 0.134, 0.136);
        TH1F *h_pi0MeanvsRun_old = new TH1F("h_pi0MeanvsRun_old", "M_{#pi0} vs Run Number;Run Number;M_{#pi0} [MeV/c^{2}]", last_run-first_run+1, first_run-0.5, last_run+0.5);
        TH1F *h_pi0SigmvsRun_old = new TH1F("h_pi0SigmvsRun_old", "#sigma_{#pi0} vs Run Number;Run Number;#sigma_{#pi0} [MeV/c^{2}]", last_run-first_run+1, first_run-0.5, last_run+0.5);
        TH1F *h_pi0MeanvsRun_new = new TH1F("h_pi0MeanvsRun_new", "M_{#pi0} vs Run Number;Run Number;M_{#pi0} [MeV/c^{2}]", last_run-first_run+1, first_run-0.5, last_run+0.5);
        TH1F *h_pi0SigmvsRun_new = new TH1F("h_pi0SigmvsRun_new", "#sigma_{#pi0} vs Run Number;Run Number;#sigma_{#pi0} [MeV/c^{2}]", last_run-first_run+1, first_run-0.5, last_run+0.5);

        TTree *t_checkCalib = new TTree("t_checkCalib", "t_checkCalib"); // a tree for checking things quicklly
        Int_t run_number;
        
        Double_t photon_px1;
        Double_t photon_py1;
        Double_t photon_pz1;
        Double_t photon_e1;
        Double_t photon_th1; // in rad
        Double_t photon_ph1; // in rad
        Double_t photon_xc1;
        Double_t photon_yc1;

        Double_t photon_px2;
        Double_t photon_py2;
        Double_t photon_pz2;
        Double_t photon_e2;
        Double_t photon_th2; // in rad
        Double_t photon_ph2; // in rad
        Double_t photon_xc2;
        Double_t photon_yc2;

        Double_t clust_x1;
        Double_t clust_y1;
        Double_t clust_ene1;
        Int_t clust_size1;

        Double_t clust_x2;
        Double_t clust_y2;
        Double_t clust_ene2;
        Int_t clust_size2;

        Double_t pi0_m;
        Double_t pi0_mm2;
        Double_t pi0_Ecal; // calculated pi0 energy
        t_checkCalib->Branch("run_number", &run_number);

        t_checkCalib->Branch("Ebeam", &Beam_energy);
        t_checkCalib->Branch("H_react_x", &H_react_x);
        t_checkCalib->Branch("H_react_y", &H_react_y);
        t_checkCalib->Branch("H_react_z", &H_react_z);
        t_checkCalib->Branch("H_gtr_px", &H_gtr_px);
        t_checkCalib->Branch("H_gtr_py", &H_gtr_py);
        t_checkCalib->Branch("H_gtr_pz", &H_gtr_pz);

        t_checkCalib->Branch("photon_px1", &photon_px1);
        t_checkCalib->Branch("photon_py1", &photon_py1);
        t_checkCalib->Branch("photon_pz1", &photon_pz1);
        t_checkCalib->Branch("photon_e1", &photon_e1);
        t_checkCalib->Branch("photon_th1", &photon_th1);
        t_checkCalib->Branch("photon_ph1", &photon_ph1);
        t_checkCalib->Branch("photon_xc1", &photon_xc1);
        t_checkCalib->Branch("photon_yc1", &photon_yc1);

        t_checkCalib->Branch("photon_px2", &photon_px2);
        t_checkCalib->Branch("photon_py2", &photon_py2);
        t_checkCalib->Branch("photon_pz2", &photon_pz2);
        t_checkCalib->Branch("photon_e2", &photon_e2);
        t_checkCalib->Branch("photon_th2", &photon_th2);
        t_checkCalib->Branch("photon_ph2", &photon_ph2);
        t_checkCalib->Branch("photon_xc2", &photon_xc2);
        t_checkCalib->Branch("photon_yc2", &photon_yc2);

        t_checkCalib->Branch("clust_x1", &clust_x1);
        t_checkCalib->Branch("clust_y1", &clust_y1);
        t_checkCalib->Branch("clust_ene1", &clust_ene1);
        t_checkCalib->Branch("clust_size1", &clust_size1); // Size of the cluster

        t_checkCalib->Branch("clust_x2", &clust_x2);
        t_checkCalib->Branch("clust_y2", &clust_y2);
        t_checkCalib->Branch("clust_ene2", &clust_ene2);
        t_checkCalib->Branch("clust_size2", &clust_size2); // Size of the cluster

        t_checkCalib->Branch("pi0_m", &pi0_m);
        t_checkCalib->Branch("pi0_mm2", &pi0_mm2);
        t_checkCalib->Branch("pi0_Ecal", &pi0_Ecal);

        // Redo clustering for each run to output their coefficients
        for(int irun = 0; irun < nRun; irun++){
            for(int iseg = 0; iseg < nSegList[irun]; iseg++){
                chain->Add(Form("%s/x%s/prodTree_%s_%d_%d.root", calibTreeDir.Data(), Kine.Data(), prefix_1.Data(), runList[irun], iseg));
                cout<<"Run:"<<runList[irun]<<"; Segment: "<<iseg<<endl;
            }

            //Read branches
            chain->SetBranchAddress("H_react_x", &H_react_x);
            chain->SetBranchAddress("H_react_y", &H_react_y);
            chain->SetBranchAddress("H_react_z", &H_react_z);
            chain->SetBranchAddress("H_gtr_px", &H_gtr_px);
            chain->SetBranchAddress("H_gtr_py", &H_gtr_py);
            chain->SetBranchAddress("H_gtr_pz", &H_gtr_pz);
            chain->SetBranchAddress("caloev", &caloev);

            // Get kinematics variables________________________________________________
            SetKineVariables(db, runList[irun]);
            
            // List of block to mask out
            caloMaskBlock = db->GetEntry_i("CALO_flag_MaskBlock", runList[irun]);
            for(int iblk = 0; iblk < nblk; iblk++){
                if(caloMaskBlock[iblk]) maskBlkList.push_back(iblk);
            }
            sort(maskBlkList.begin(), maskBlkList.end()); // sort the block numbers
            auto last = unique(maskBlkList.begin(), maskBlkList.end()); // move the duplicate numbers to the end
            maskBlkList.erase(last, maskBlkList.end());  // remove the duplicate numbers

            // Beam electron and target 4-momenta for missing mass calculation
            TLorentzVector beam(0, 0, Beam_energy, Beam_energy);
            TLorentzVector p0(0, 0, 0, M_targ);
             
            // initialize variables for saving tree
            run_number = runList[irun];
            photon_px1 = -999;
            photon_py1 = -999;
            photon_pz1 = -999;
            photon_e1 = -999;
            photon_th1 = -999;
            photon_ph1 = -999;
            photon_xc1 = -999;
            photon_yc1 = -999;

            photon_px2 = -999;
            photon_py2 = -999;
            photon_pz2 = -999;
            photon_e2 = -999;
            photon_th2 = -999;
            photon_ph2 = -999;
            photon_xc2 = -999;
            photon_yc2 = -999;

            clust_x1 = -999;
            clust_y1 = -999;
            clust_ene1 = -999;
            clust_size1 = -999;
            
            clust_x2 = -999;
            clust_y2 = -999;
            clust_ene2 = -999;
            clust_size2 = -999;
            
            pi0_m = -999;
            pi0_mm2 = -999;
            pi0_Ecal = -999;

            // Event loop: check the mass peak with old coefficients
            nevt = chain->GetEntries();
            for(Int_t i = 0; i < nevt; i++){ // Event loop: check the mass peak with old coefficients
                Int_t ObjectNumber = TProcessID::GetObjectCount();
                chain->GetEntry(i);
                if (i % 1000 == 0) cout << i << "/" << nevt << endl;

                //Redo the cluster to update the pi0 mass 
                for(int iblk = 0; iblk < nblk; iblk++) hasBlock[iblk] = false;

                Int_t nCaloBlock = caloev->GetNbBlocks();
                for(int iblk = 0; iblk < nCaloBlock; iblk++){
                    TCaloBlock *block = caloev->GetBlock(iblk);
                    Int_t nb = block->GetBlockNumber(); // Block number in simulation numbering scheme
                    Float_t energy = block->GetEnergy(0)*pi0_coef[bnConv_OldToNew(nb)]/coefElas[bnConv_OldToNew(nb)]; // Replace the coefficient

                    block->Erase("");
                    block->AddPulse(energy, 0);
                    hasBlock[nb] = true;
                }
                
                for(Int_t iblk = 0; iblk < nblk; iblk++){
                    if(hasBlock[iblk] == false) caloev->AddBlock(iblk);
                }

                // caloev->TriggerSim(0.2);
                caloev->DoClustering(-3,3); // Time window (-1.5,1.5) [ns]
                Int_t nclus = caloev->GetNbClusters();
                // cout<<nclus<<endl;

                for (Int_t iclus = 0; iclus < nclus; iclus++){
                    caloev->GetCluster(iclus)->Analyze(1, -3, 3);
                }

                if(nclus != 2) continue;

                if(!(TMath::Abs(H_react_z) < 15)) continue;
                ev->SetVertex(H_react_x, H_react_y, H_react_z);
                ev->GetGeometry()->SetCaloTheta(-1*NPS_angle);
                ev->GetGeometry()->SetCaloDist(NPS_dist);
                ev->SetCaloEvent(caloev);

                // Check seed blocks and apply acceptance cut
                Int_t seed1 = GetSeedBlocks(caloev->GetCluster(0)->GetX(), caloev->GetCluster(0)->GetY());
                Int_t seed2 = GetSeedBlocks(caloev->GetCluster(1)->GetX(), caloev->GetCluster(1)->GetY());
                if(!(PassAccCut(maskBlkList, NPS_dist, seed1) && PassAccCut(maskBlkList, NPS_dist, seed2))) continue;

                // Get reconstructed photons
                TLorentzVector photon1 = ev->GetPhoton(0, 7, 0);
                TLorentzVector photon2 = ev->GetPhoton(1, 7, 0);
                Double_t m_pi0_recons = (photon1 + photon2).M();
                
                TLorentzVector kpvec;
                kpvec.SetXYZM(H_gtr_px, H_gtr_py, H_gtr_pz, 0.000511);
                Double_t mm2_pi0_recons = (beam + p0 - kpvec - photon1 - photon2).Mag2();

                if(photon1.E() > photEcut && photon2.E() > photEcut) {
                    h_pi0M_final->Fill(m_pi0_recons);
                    h_pi0Mx2_final->Fill(mm2_pi0_recons);
                }
                ev->Reset();
                TProcessID::SetObjectCount(ObjectNumber);
            } // loop of events for old mass plots 
            
            f_fit_final->SetParLimits(1, 0.1, 0.15);
            f_fit_final->SetParLimits(2, 0.00001, 0.03);
            Double_t meanfit = h_pi0M_final->GetBinCenter(h_pi0M_final->GetMaximumBin());
            Double_t minfit = meanfit - 0.035;
            Double_t maxfit = meanfit + 0.035;
            h_pi0M_final->Fit(f_fit_final, "R", "", minfit, maxfit);

            N = f_fit_final->GetParameter(0)/h_pi0M_final->GetBinWidth(1);
            N_err = f_fit_final->GetParError(0)/h_pi0M_final->GetBinWidth(1);
            mean = f_fit_final->GetParameter(1);
            mean_err = f_fit_final->GetParError(1);
            sigm = f_fit_final->GetParameter(2);
            sigm_err = f_fit_final->GetParError(2);

            Int_t ibin = runList[irun]-first_run+1;
            h_pi0MeanvsRun_old->SetBinContent(ibin, mean*1000);
            h_pi0MeanvsRun_old->SetBinError(ibin, mean_err*1000);
            h_pi0SigmvsRun_old->SetBinContent(ibin, sigm*1000);
            h_pi0SigmvsRun_old->SetBinError(ibin, sigm_err*1000);

            outfile->cd();
            h_pi0M_final->Write(Form("h_pi0M_old_%d", runList[irun]));
            h_pi0M_final->Reset();
            h_pi0Mx2_final->Write(Form("h_pi0Mx2_old_%d", runList[irun]));
            h_pi0Mx2_final->Reset();

            // Update new coefficients for this run
            Double_t pi0_coef_new[nblk];
            for(int iblk = 0; iblk < nblk; iblk++) pi0_coef_new[iblk] = 0;
            for(int iblk = 0; iblk < nblk; iblk++) pi0_coef_new[iblk] = pi0_coef[iblk]*(m_pi0 / mean); // in NPS numbering scheme

            // output to text file
            system(Form("rm -f Result/%s/coef_pi0Calib_%d.txt", filename.Data(), runList[irun]));
            ofstream fcoef_out(Form("Result/%s/coef_pi0Calib_%d.txt", filename.Data(), runList[irun]));
            for(int iblk = 0; iblk < nblk; iblk++) fcoef_out<<pi0_coef_new[iblk]<<endl;
            fcoef_out.close();

            // output to text file for nps replay
            system(Form("rm -f Result/%s/coef_pi0Calib_%d_npsreplay.txt", filename.Data(), runList[irun]));
            ofstream fout_coef_hcana(Form("Result/%s/coef_pi0Calib_%d_npsreplay.txt", filename.Data(), runList[irun]));
            fout_coef_hcana<<"nps_cal_arr_gain_cor = ";
            for(int iblk = 0; iblk < nblk; iblk++){
                if((iblk+1) % 30 == 0) fout_coef_hcana<<pi0_coef_new[iblk]<<", \n           ";
                else fout_coef_hcana<<pi0_coef_new[iblk]<<", ";
            }

            // Event loop: check the mass peak with updated coefficients
            for(Int_t i = 0; i < nevt; i++){ // Event loop: check the mass peak with updated coefficients
                Int_t ObjectNumber = TProcessID::GetObjectCount();
                chain->GetEntry(i);
                if (i % 1000 == 0) cout << i << "/" << nevt << endl;
                
                //Redo the cluster to update the pi0 mass 
                for(int iblk = 0; iblk < nblk; iblk++) hasBlock[iblk] = false;

                Int_t nCaloBlock = caloev->GetNbBlocks();
                for(int iblk = 0; iblk < nCaloBlock; iblk++){
                    TCaloBlock *block = caloev->GetBlock(iblk);
                    Int_t nb = block->GetBlockNumber(); // Block number in simulation numbering scheme
                    Float_t energy = block->GetEnergy(0)*pi0_coef_new[bnConv_OldToNew(nb)]/coefElas[bnConv_OldToNew(nb)]; // Replace the coefficient

                    block->Erase("");
                    block->AddPulse(energy, 0);
                    hasBlock[nb] = true;
                }
                
                for(Int_t iblk = 0; iblk < nblk; iblk++){
                    if(hasBlock[iblk] == false) caloev->AddBlock(iblk);
                }

                // caloev->TriggerSim(0.2);
                caloev->DoClustering(-3,3); // Time window (-1.5,1.5) [ns]
                Int_t nclus = caloev->GetNbClusters();
                // cout<<nclus<<endl;

                for (Int_t iclus = 0; iclus < nclus; iclus++){
                    caloev->GetCluster(iclus)->Analyze(1, -3, 3);
                }

                if(nclus != 2) continue;

                if(!(TMath::Abs(H_react_z) < 15)) continue;
                ev->SetVertex(H_react_x, H_react_y, H_react_z);
                ev->GetGeometry()->SetCaloTheta(-1*NPS_angle);
                ev->GetGeometry()->SetCaloDist(NPS_dist);
                ev->SetCaloEvent(caloev);

                // Check seed blocks and apply acceptance cut
                Int_t seed1 = GetSeedBlocks(caloev->GetCluster(0)->GetX(), caloev->GetCluster(0)->GetY());
                Int_t seed2 = GetSeedBlocks(caloev->GetCluster(1)->GetX(), caloev->GetCluster(1)->GetY());
                if(!(PassAccCut(maskBlkList, NPS_dist, seed1) && PassAccCut(maskBlkList, NPS_dist, seed2))) continue;

                // Get the cluster position and reconstructed photon
                clust_x1   = caloev->GetCluster(0)->GetX();
                clust_y1   = caloev->GetCluster(0)->GetY();
                clust_ene1 = caloev->GetCluster(0)->GetE();
                clust_size1   = caloev->GetCluster(0)->GetClusSize();

                clust_x2   = caloev->GetCluster(1)->GetX();
                clust_y2   = caloev->GetCluster(1)->GetY();
                clust_ene2 = caloev->GetCluster(1)->GetE();
                clust_size2   = caloev->GetCluster(1)->GetClusSize();

                // photon informations
                TLorentzVector photon1 = ev->GetPhoton(0, 7, 0);
                TLorentzVector photon2 = ev->GetPhoton(1, 7, 0);

                photon_px1 = photon1.Px();
                photon_py1 = photon1.Py();
                photon_pz1 = photon1.Pz();
                photon_e1 = photon1.E();
                photon_th1 = photon1.Theta();
                photon_ph1 = photon1.Phi();

                photon_px2 = photon2.Px();
                photon_py2 = photon2.Py();
                photon_pz2 = photon2.Pz();
                photon_e2 = photon2.E();
                photon_th2 = photon2.Theta();
                photon_ph2 = photon2.Phi();

                // Impact position of reconstructed photons on NPS surface
                Double_t beta = TMath::ATan(photon1.Px() / photon1.Pz());
                photon_xc1 = NPS_dist*TMath::Tan(beta-NPS_angle)-H_react_z*TMath::Sin(beta)/TMath::Cos(beta-NPS_angle);
                photon_yc1 = (TMath::Sqrt(NPS_dist*NPS_dist+photon_xc1*photon_xc1)-H_react_z*(TMath::Cos(beta)+TMath::Sin(beta)*TMath::Tan(beta-NPS_angle))) * photon1.Py();
                photon_yc1 = photon_yc1/TMath::Sqrt(photon1.Px()*photon1.Px()+photon1.Pz()*photon1.Pz());

                beta = TMath::ATan(photon2.Px() / photon2.Pz());
                photon_xc2 = NPS_dist*TMath::Tan(beta-NPS_angle)-H_react_z*TMath::Sin(beta)/TMath::Cos(beta-NPS_angle);
                photon_yc2 = (TMath::Sqrt(NPS_dist*NPS_dist+photon_xc2*photon_xc2)-H_react_z*(TMath::Cos(beta)+TMath::Sin(beta)*TMath::Tan(beta-NPS_angle)))*photon2.Py();
                photon_yc2 = photon_yc2/TMath::Sqrt(photon2.Px()*photon2.Px()+photon2.Pz()*photon2.Pz());

                Double_t H_gtr_e = sqrt(H_gtr_px*H_gtr_px+H_gtr_py*H_gtr_py+H_gtr_pz*H_gtr_pz+0.000511*0.000511);
                TLorentzVector Lk = beam;
                TLorentzVector Lp = p0;
                TLorentzVector Lkp(H_gtr_px, H_gtr_py, H_gtr_pz, H_gtr_e);

                // pi0 invariant mass and missing mass
                pi0_m = (photon1 + photon2).M();
                pi0_mm2 = (beam + p0 - Lkp - photon1 - photon2).Mag2();

                // Expected pi0 energy
                TVector3 Vk, Vkp, Vq1, Vq2;
                Vk = Lk.Vect();
                Vkp = Lkp.Vect();
                Vq1 = photon1.Vect();
                Vq2 = photon2.Vect();

                Double_t E_kkpN, P_vq, costheta;
                E_kkpN = Lk.E()-Lkp.E()+M_targ; // energy of E_k-E_kp+M_targ
                P_vq = (Vk-Vkp).Mag(); // magnitude of the virtual photon momentum
                costheta = TMath::Cos((Vk-Vkp).Angle(Vq1+Vq2)); // angle between virtual photon and pi0

                Double_t a, b, c;
                a = 4*E_kkpN*E_kkpN-4*P_vq*P_vq*costheta*costheta;
                b = 4*E_kkpN*(M_targ*M_targ-(Lk-Lkp+Lp).M2()-m_pi0*m_pi0);
                c = 4*m_pi0*m_pi0*P_vq*P_vq*costheta*costheta+(M_targ*M_targ-(Lk-Lkp+Lp).M2()-m_pi0*m_pi0)*(M_targ*M_targ-(Lk-Lkp+Lp).M2()-m_pi0*m_pi0);
                pi0_Ecal = (-1*b+sqrt(b*b-4*a*c))/(2*a); // predicted pi0 energy

                if(photon1.E() > photEcut && photon2.E() > photEcut) {
                    h_pi0M_final->Fill(pi0_m);
                    h_pi0Mx2_final->Fill(pi0_mm2);
                }

                t_checkCalib->Fill();
                ev->Reset();
                TProcessID::SetObjectCount(ObjectNumber);
            } // loop of events updated mass plot

            f_fit_final->SetParLimits(1, 0.1, 0.15);
            f_fit_final->SetParLimits(2, 0.00001, 0.03);
            meanfit = h_pi0M_final->GetBinCenter(h_pi0M_final->GetMaximumBin());
            minfit = meanfit - 0.035;
            maxfit = meanfit + 0.035;
            h_pi0M_final->Fit(f_fit_final, "R", "", minfit, maxfit);

            N = f_fit_final->GetParameter(0)/h_pi0M_final->GetBinWidth(1);
            N_err = f_fit_final->GetParError(0)/h_pi0M_final->GetBinWidth(1);
            mean = f_fit_final->GetParameter(1);
            mean_err = f_fit_final->GetParError(1);
            sigm = f_fit_final->GetParameter(2);
            sigm_err = f_fit_final->GetParError(2);

            h_pi0MeanvsRun_new->SetBinContent(ibin, mean*1000);
            h_pi0MeanvsRun_new->SetBinError(ibin, mean_err*1000);
            h_pi0SigmvsRun_new->SetBinContent(ibin, sigm*1000);
            h_pi0SigmvsRun_new->SetBinError(ibin, sigm_err*1000);

            outfile->cd();
            h_pi0M_final->Write(Form("h_pi0M_new_%d", runList[irun]));
            h_pi0M_final->Reset();
            h_pi0Mx2_final->Write(Form("h_pi0Mx2_new_%d", runList[irun]));
            h_pi0Mx2_final->Reset();

            chain->Reset();
        } // loop of runs

        outfile->cd();
        t_checkCalib->Write();
        h_coEff->Write();
        hh_coEff->Write();
        h_pi0MeanvsRun_old->Write();
        h_pi0MeanvsRun_new->Write();
        h_pi0SigmvsRun_old->Write();
        h_pi0SigmvsRun_new->Write();
        outfile->Close();
    } // if iIter=9
}

bool PassAccCut(const vector<int> &maskBlkList, const double &NPS_dist, int sim_blk_number)
{
    bool pass = true; // assume all pass for initialization

    // Column and row numbers for the cut on seed blocks
    // these col and row numbers are in sim. numbering scheme
    Int_t sim_irow = sim_blk_number%nrow;
    Int_t sim_icol = (sim_blk_number-sim_irow)/nrow;

    // shell not pass if at the edge columns/rows====================
    if(!(0 < sim_irow && sim_irow < 35 && 0 < sim_icol && sim_icol < 29)){
        pass = false; 
        return pass;
    }

    // shell not pass if next to a bad/off block=============================
    Int_t nps_blk_number = bnConv_OldToNew(sim_blk_number); // Convert block number, column, row to NPS number scheme
    Int_t icol = nps_blk_number%30;
    Int_t irow = (nps_blk_number-icol)/30;
    // Check if next to masked block
    for(int i = 0; i < maskBlkList.size(); i++){
        Int_t icol_maskBlk = maskBlkList[i]%30;
        Int_t irow_maskBlk = (maskBlkList[i]-icol)/30;
        if(abs(icol - icol_maskBlk) <= 1
        && abs(irow - irow_maskBlk) <= 1){
            pass = false; // shell not pass if next to a masked block
            return pass;
        }
    }

    // shell not pass if under the shadow of magnet (for better statistics)===========
    Int_t nCol_inMgnt = 0; // number of columns under the shadow of magnet
    if(290 < NPS_dist && NPS_dist < 310) nCol_inMgnt = 5; // col0-4 when NPS@300cm
    else if(340 < NPS_dist && NPS_dist < 360) nCol_inMgnt = 4; // col0-3 when NPS@350cm
    else if(390 < NPS_dist && NPS_dist < 410) nCol_inMgnt = 2; // col0-1 when NPS@400cm
    else nCol_inMgnt = 0; // no blocked when NPS is far from target
    if(!(abs(sim_icol-29) >= nCol_inMgnt)){
        pass = false; // shell not pass if at the edge columns/rows
        return pass;
    }

    return pass;
}

void SetKineVariables(TDVCSDB *db, int run_number){
    Beam_energy = *db->GetEntry_d("BEAM_param_Energy", run_number);
    HMS_mom = *db->GetEntry_d("SIMU_param_HMSmomentum", run_number);
    HMS_angle = *db->GetEntry_d("SIMU_param_HMSangle", run_number);
    Target_amu = *db->GetEntry_d("TARGET_param_Amu", run_number);

    // The angle and distance of calorimeter
    NPS_dist = *db->GetEntry_d("CALO_geom_Dist", run_number);
    NPS_angle = *db->GetEntry_d("CALO_geom_Yaw", run_number); // unit in rad already in DB

    // Elastic coefficients
    coefElas = db->GetEntry_d("CALO_calib_ElasCoef", run_number);
}

// void GetSeedBlocks(TCaloEvent *caloev, Int_t &seed1, Int_t &seed2){
//     Double_t maxene;
//     for(int iclus = 0; iclus < 2; iclus++){
//         maxene = -1;
//         if(iclus == 0) seed1 = -1;
//         else seed2 = -1;
//         for (Int_t k = 0; k < caloev->GetCluster(iclus)->GetClusSize(); k++){
//             TCaloBlock *clusblock = caloev->GetCluster(iclus)->GetBlock(k);
//             Double_t blockene = clusblock->GetBlockEnergy();
//             if (blockene > maxene){
//                 maxene = blockene;
//                 if(iclus == 0) seed1 = clusblock->GetBlockNumber();
//                 else seed2 = clusblock->GetBlockNumber();
//             }
//         }
//     }
//     // Check if the seeds are reasonable
//     if(seed1 <  0 || seed2 < 0 || seed1 > 1079 || seed2 > 1079){
//         cout<<"ERROR: Got a wrong seed block!!!"<<endl;
//         return;
//     }
// }

int GetSeedBlocks(Double_t xclus, Double_t yclus){
    // determine the column and row numbers of the cluster
    Double_t step  = xnut/15; //  = ynut/18
    Double_t xclus_raw = xnut - xclus;
    Double_t yclus_raw = ynut + yclus;
    Int_t sim_icol = std::round((xclus_raw - step/2) / step); // column in simulation numbering scheme
    Int_t sim_irow = std::round((yclus_raw - step/2) / step); // row in simulation numbering scheme
    Int_t sim_blk_number = sim_icol*nrow + sim_irow; // block number in simulation numbering scheme
    // cout<<step<<" "<<xclus<<" "<<yclus<<" "<<sim_icol<<" "<<sim_irow<<" "<<sim_blk_number<<endl;

    // Check if the seeds are reasonable
    if(sim_blk_number <  0 || sim_blk_number > 1079){
        cout<<"ERROR: Got a wrong seed block!!!"<<endl;
        sim_blk_number = -999;
    }
    return sim_blk_number;
}