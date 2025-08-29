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
#include <unistd.h>  // for getpid()

// Number of columns used for calibration
const Int_t ncol_inUse = 23;
const Int_t nblk_inUse = ncol_inUse*nrow;

// Check if the block is in the range for calibration (with sim. numbering scheme)
bool passAccCut(int sim_blk_number);

void PrintMemoryUsage(const std::string& tag = "") {
    pid_t pid = getpid();
    std::ifstream statm("/proc/self/statm");

    if (!statm.is_open()) {
        std::cerr << "Failed to open /proc/self/statm" << std::endl;
        return;
    }

    long size = 0, resident = 0, shared = 0;
    statm >> size >> resident >> shared; // Memory usage, unit in pages

    long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // page size = 4096 bytes = 4 kb
    long rss = resident * page_size_kb;
    long vsz = size * page_size_kb;

    if(rss < 1000) std::cout << "[Memory] " << tag << "  RSS = " << rss << " KB" << " | VSZ = " << vsz << " KB" << std::endl;
    else if(1000 < rss && rss < 1000000) std::cout << "[Memory] " << tag << "  RSS = " << rss/1000 << " MB" << " | VSZ = " << vsz/1000 << " MB" << std::endl;
    if(rss > 1000000) std::cout << "[Memory] " << tag << "  RSS = " << rss/1000000 << " GB" << " | VSZ = " << vsz/1000000 << " GB" << std::endl;
}

/*double crystalball_function(double *x, double *p)
{
    // see math/mathcore/src/PdfFuncMathCore.cxx in ROOT 6.x
    double N = p[0];
    double mean = p[1];
    double sigma = p[2];
    double alpha = p[3];
    double n = p[4];
    if ((!x) || (!p)) return 0.; // just a precaution
    
    // evaluate the crystal ball function
    if (sigma < 0.)     return 0.;
    double z = (x[0] - mean)/sigma;
    if (alpha < 0) z = -z; 
    double abs_alpha = std::abs(alpha);
    // double C = n/abs_alpha * 1./(n-1.) * std::exp(-alpha*alpha/2.);
    // double D = std::sqrt(M_PI/2.)*(1.+ROOT::Math::erf(abs_alpha/std::sqrt(2.)));
    // double N = 1./(sigma*(C+D));
    if (z  > - abs_alpha)
      return N * std::exp(- 0.5 * z * z);
    else {
      //double A = std::pow(n/abs_alpha,n) * std::exp(-0.5*abs_alpha*abs_alpha);
      double nDivAlpha = n/abs_alpha;
      double AA =  std::exp(-0.5*abs_alpha*abs_alpha);
      double B = nDivAlpha -abs_alpha;
      double arg = nDivAlpha/(B-z);
      return N * AA * std::pow(arg,n);
    }
}

double polynomial_function(double *x, double *p)
{
    return p[0]+p[1]*x[0];
}

double crystalball_poly(double *x, double *p)
{
    return crystalball_function(x, p) + polynomial_function(x, &p[5]);
}*/

void pi0Calib(TString Kine, int TargetFlag, int iIter)
{   
    gSystem->Load("/group/nps/hhuang/software/NPS_SOFT/libDVCS.so");
    TDVCSDB *db = new TDVCSDB("dvcs", "clrlpc", 3306, "hhuang", "");

    TString Tar;
    if(TargetFlag == 0) Tar = "LH2";
    else if(TargetFlag == 1) Tar = "LD2";
    else if(TargetFlag == -1) Tar = "LH2_LD2";
    else{
        cout<<"ERROR: Unknown target!!!"<<endl;
        return;
    }

    // Make folder for the new calibration results
    TString filename = Form("%s_%s", Kine.Data(), Tar.Data());
    if(iIter == 1) system(Form("mkdir Result/%s_pass2_v3", filename.Data()));

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

    // for(int irun = 0; irun < nRun; irun++){
    //     cout<<"Run: "<<runList[irun]<<"; nSeg: "<<nSegList[irun]<<endl;
    // }

    // Add the rootfiles to TChain
    TChain *chain = new TChain("t_prod");
    for(int irun = 0; irun < nRun; irun++){
        for(int iseg = 0; iseg < nSegList[irun]; iseg++){
            chain->Add(Form("/group/nps/hhuang/analysis/DVCS_NPS2023/DVCS_analysis/pi0Calib_wf/calibTree/%s/prodTree_pass2_v3_%d_%d.root", Kine.Data(), runList[irun], iseg));
            // cout<<"Run:"<<runList[irun]<<"; Segment: "<<iseg<<endl;
        }
    }
    fRunList.close();

    // Get branches in the rootfiles
    Double_t hms_Vx;
    Double_t hms_Vy;
    Double_t hms_Vz;

    TCaloEvent *caloev = new TCaloEvent();

    chain->SetBranchAddress("hms_Vx", &hms_Vx);
    chain->SetBranchAddress("hms_Vy", &hms_Vy);
    chain->SetBranchAddress("hms_Vz", &hms_Vz);
    chain->SetBranchAddress("caloev", &caloev);

    // Get kinematics variables________________________________________________
    Double_t Beam_energy = *db->GetEntry_d("BEAM_param_Energy", first_run);
    Double_t HMS_mom = *db->GetEntry_d("SIMU_param_HMSmomentum", first_run);
    Double_t HMS_angle = *db->GetEntry_d("SIMU_param_HMSangle", first_run);
    Double_t Target_amu = *db->GetEntry_d("TARGET_param_Amu", first_run);

    // The angle and distance of calorimeter
    Double_t NPS_dist = *db->GetEntry_d("CALO_geom_Dist", first_run);
    Double_t NPS_angle = *db->GetEntry_d("CALO_geom_Yaw", first_run); // Rad already in DB
    Int_t *caloMaskBlock = new Int_t[1080]; // Get the mask block information in NPS numbering Scheme
    caloMaskBlock = db->GetEntry_i("CALO_flag_MaskBlock", first_run);

    // Elastic coefficients
    Double_t *coefElas = new Double_t[1080]; // Get the elastic coefficients in NPS numbering Scheme
    coefElas = db->GetEntry_d("CALO_calib_ElasCoef", first_run);

    // for (Int_t i = 0; i < 1080; i++) cout<<coef[i]<<endl;

    // Print out some basic information for recording
    system(Form("rm -f Result/%s_pass2_v3/calib_Info.txt", filename.Data())); // avoid appending old file
    ofstream info_Stream(Form("Result/%s_pass2_v3/calib_Info.txt", filename.Data()));

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
    info_Stream<<"Run list: calibList/"<<Kine.Data()<<"/"<<Tar.Data()<<".txt"<<endl;
    for(int irun = 0; irun < nRun; irun++) info_Stream<<"Run "<<runList[irun]<<endl;
    info_Stream<<" "<<endl;

    // Initial settings for calibration
    TLorentzVector beam(0, 0, Beam_energy, Beam_energy);
    TLorentzVector p0(0, 0, 0, m_p);

    TDVCSEvent *ev = new TDVCSEvent();
    Float_t alpha = NPS_angle, d = NPS_dist;
    Bool_t hasBlock[1080];

    // matrix and vectors used to compute epsilon
    Float_t B = 0;
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

    // Ouput files for recording
    system(Form("rm -f Result/%s_pass2_v3/log_Iteration_%d.txt", filename.Data(), iIter)); // avoid appending old file
    ofstream log_Stream(Form("Result/%s_pass2_v3/log_Iteration_%d.txt", filename.Data(), iIter));
    
    log_Stream<<"++++++++++ This is iteration "<<iIter<<" of calibration ++++++++++"<<endl;

    Int_t nevt;
    cout<<chain->GetEntries()<<endl;
    if(chain->GetEntries() >= 200000) nevt = 200000;
    else nevt = chain->GetEntries();

    if(iIter == 1){ // initialize the correction factor and generate the mass distribution from clusters
        for(Int_t i = 0; i < nblk; i++) corr_pi0[i] = 1; // initialize the correction factor

        // Output histogram of invariant mass distribution
        TFile *output_iter0 = new TFile(Form("Result/%s_pass2_v3/f_Mgg_%d.root", filename.Data(), iIter-1), "recreate");
        TH1F *h_pi0M_iter0 = new TH1F(Form("h_pi0M_%d", iIter-1), "M_{#gamma#gamma} (E_{#gamma1} > 1.4 GeV && E_{#gamma2} > 1.4 GeV);M_{#gamma#gamma} [GeV];Counts", 150, 0.05, 0.2);
        h_pi0M_iter0->Sumw2();
        // TF1 *f_fit_iter0 = new TF1(Form("f_fit_%d", iIter-1), "gausn(0)+pol1(3)", 0.08, 0.2);
        // TF1 *f_fit_iter0 = new TF1(Form("f_fit_%d", iIter-1), crystalball_poly, 0.08, 0.2, 7);

        for(Int_t i = 0; i < nevt; i++){ // Event loop: get the mass distribution from clusters
            chain->GetEntry(i);
            // if (i % 1000 == 0) cout << i << "/" << nevt << endl;
            if(i % 10000 == 0) PrintMemoryUsage(Form("After processing %d events", i));

            ev->SetVertex(hms_Vx, hms_Vy, hms_Vz);
            ev->GetGeometry()->SetCaloTheta(-1*NPS_angle);
            ev->GetGeometry()->SetCaloDist(NPS_dist);
            ev->SetCaloEvent(caloev);

            if(caloev->GetNbClusters() != 2) continue; // Only look at two clusters events

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
            m_pi0_recons = (photon1 + photon2).M();
            if(photon1.E() > 1.4 && photon2.E() > 1.4) h_pi0M_iter0->Fill(m_pi0_recons);

            ev->Reset();
        }// event loop

        Double_t meanfit_iter0 = h_pi0M_iter0->GetBinCenter(h_pi0M_iter0->GetMaximumBin());
        Double_t minfit_iter0 = meanfit_iter0 - 0.01;
        Double_t maxfit_iter0 = meanfit_iter0 + 0.01;
        TF1 *f_fit_iter0 = new TF1(Form("f_fit_%d", iIter-1), "gausn(0)", minfit_iter0, maxfit_iter0);
        // f_fit_iter0->SetParLimits(0, 100, 100000);
        f_fit_iter0->SetParLimits(1, 0.115, 0.16);
        f_fit_iter0->SetParLimits(2, 0.00001, 0.03);
        // f_fit_iter0->SetParameters(10000, 0.13, 0.001, 1.2, 1.0); // N, mean, sigma, alpha, n
        // f_fit_iter0->SetParLimits(2, 0.00001, 0.03);
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
     
    else{ // Get the correction factor and mass distrubution from the previous iteration when iIter > 1
        ifstream fcorr_old(Form("Result/%s_pass2_v3/corr_pi0.txt", filename.Data()));
        for(int iblk = 0; iblk < 1080; iblk++) fcorr_old>>corr_pi0[iblk];

        // Get the mass distribution with the histogram from previous iteration
        TFile *input = TFile::Open(Form("Result/%s_pass2_v3/f_Mgg_%d.root", filename.Data(), iIter-1));
        TH1F *h_pi0M_last = (TH1F*)input->Get(Form("h_pi0M_%d", iIter-1));
        // TF1 *f_fit_last = new TF1("f_fit_last", "gausn(0)+pol1(3)", 0.08, 0.2);
        // TF1 *f_fit_last = new TF1("f_fit_last", crystalball_poly, 0.08, 0.2, 7);

        Double_t meanfit_last = h_pi0M_last->GetBinCenter(h_pi0M_last->GetMaximumBin());
        Double_t minfit_last = meanfit_last - 0.01;
        Double_t maxfit_last = meanfit_last + 0.01;
        TF1 *f_fit_last = new TF1(Form("f_fit_%d", iIter-1), "gausn(0)", minfit_last, maxfit_last);
        // f_fit_last->SetParLimits(0, 100, 100000);
        f_fit_last->SetParLimits(1, 0.12, 0.15);
        f_fit_last->SetParLimits(2, 0.00001, 0.03);
        // f_fit_last->SetParameters(10000, 0.13, 0.001, 1.2, 1.0); // N, mean, sigma, alpha, n
        // f_fit_last->SetParLimits(2, 0.00001, 0.03);
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
    }

    log_Stream<<"+++++++++++++++++++++++++"<<endl;
    log_Stream<<" "<<endl;
    log_Stream<<"Start the pi0 calibration......"<<endl;
    
    for (Int_t i = 0; i < nevt; i++){ // Event loop: pi0 calibration
        chain->GetEntry(i);
        // if (i % 1000 == 0) cout << i << "/" << nevt << endl;
        if(i % 10000 == 0) PrintMemoryUsage(Form("After processing %d events", i));

        //Redo the cluster to update the pi0 mass
        for(int iblk = 0; iblk < 1080; iblk++) hasBlock[iblk] = false;

        Int_t nCaloBlock = caloev->GetNbBlocks();
        for(int iblk = 0; iblk < nCaloBlock; iblk++){
            TCaloBlock *block = caloev->GetBlock(iblk);
            Int_t nb = block->GetBlockNumber();
            Float_t energy = block->GetEnergy(0)*corr_pi0[nb];

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
        if(!(caloev->GetCluster(0)->GetEnergy() > 1.4 && caloev->GetCluster(1)->GetEnergy() > 1.4)) continue;

        ev->SetVertex(hms_Vx, hms_Vy, hms_Vz);
        ev->GetGeometry()->SetCaloTheta(-1*NPS_angle);
        ev->GetGeometry()->SetCaloDist(NPS_dist);
        ev->SetCaloEvent(caloev);

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
        
        // Calibration
        B += m_pi0_recons*m_pi0_recons - m_pi0*m_pi0;

        // compute C_mat matrix and vectors D and L
        for (Int_t j = 0; j < nblk; j++){
            dmdepsilon[j] = m_pi0_recons * m_pi0_recons * block_energy_ratio[j];
        }

        for(Int_t j = 0; j < nblk; j++){
            Int_t jbn_sim = j;
            if (dmdepsilon[jbn_sim] != 0){
                D[j] -= (m_pi0_recons * m_pi0_recons - m_pi0 * m_pi0) * dmdepsilon[jbn_sim];
                L[j] += dmdepsilon[jbn_sim];

                for (Int_t k = 0; k < nblk; k++){
                    Int_t kbn_sim = k;
                    C_mat(j, k) += dmdepsilon[jbn_sim] * dmdepsilon[kbn_sim];
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
    } // Event loop

    // Calculate calibration coefficients
    Int_t count_zero = 0;
    for(Int_t j = 0; j < nblk; j++){
        for(Int_t k = 0; k < nblk; k++){
            if(j != k) continue;
            if(C_mat(j, k) == 0){
                cout<<bnConv_OldToNew(j)<<endl;
                count_zero += 1;
            }
        }
    }

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
    
    // computation of epsilons
    for (Int_t i = 0; i < nblk; i++){
        for (Int_t j = 0; j < nblk; j++){
            LCD += C_mat_inv(i, j) * D[j] * L[i];
            LCL += C_mat_inv(i, j) * L[j] * L[i];
        }
    }

    lambda = (B + LCD) / (LCL);

    for (Int_t i = 0; i < nblk; i++){
        for (Int_t j = 0; j < nblk; j++){
            epsilon[i] += C_mat_inv(i, j) * (D[j] - lambda * L[j]);
        }
    }

    // update correction factors with epsilon
    Double_t mean_pi0 = 0.5*(pi0Mass_min+pi0Mass_max);
    for (Int_t i = 0; i < nblk; i++) corr_pi0[i] *= (1 + epsilon[i]);

    // update correction coefficients for those outside acceptance with the average of the others
    Double_t sum_inAcc = 0;
    Double_t count_inAcc = 0;
    for (Int_t i = 0; i < 1080; i++){
        if(passAccCut(i)){
            sum_inAcc+=coefElas[bnConv_OldToNew(i)]*corr_pi0[i];
            count_inAcc+=1.;
        }
        // cout<<corr_pi0[i]<<" : "<<sum_inAcc<<endl;
    }

    for (Int_t i = 0; i < 1080; i++){
        if(!passAccCut(i)) corr_pi0[i] = (sum_inAcc/count_inAcc)/coefElas[bnConv_OldToNew(i)];
    }

    // save corr_pi0 coeffs
    system(Form("rm -f Result/%s_pass2_v3/corr_pi0.txt", filename.Data())); // avoid appending old file
    ofstream corr_pi0_stream(Form("Result/%s_pass2_v3/corr_pi0.txt", filename.Data()));
    for (Int_t i = 0; i < 1080; i++){
        corr_pi0_stream<<corr_pi0[i]<<endl;
        // cout << corr_pi0[i] << endl;
    }
    corr_pi0_stream.close();

    log_Stream<<"The pi0 calibration of iteration "<<iIter<<" finished"<<endl;
    log_Stream<<"Update the invariant mass distribution with last iteration......"<<endl;
    
    // Event loop: check the new mass position of pi0
    TFile *output = new TFile(Form("Result/%s_pass2_v3/f_Mgg_%d.root", filename.Data(), iIter), "recreate");
    TH1F *h_pi0M = new TH1F(Form("h_pi0M_%d", iIter), "M_{#gamma#gamma} (E_{#gamma1} > 1.4 GeV && E_{#gamma2} > 1.4 GeV);M_{#gamma#gamma} [GeV];Counts", 150, 0.05, 0.2);
    h_pi0M->Sumw2();
    // TF1 *f_fit = new TF1(Form("f_fit_%d", iIter), "gausn(0)+pol1(3)", 0.08, 0.2);
    // TF1 *f_fit = new TF1(Form("f_fit_%d", iIter), crystalball_poly, 0.08, 0.2, 7);

    for (Int_t i = 0; i < nevt; i++){ // Event loop: update the invariant mass distribution after the last iteration
        chain->GetEntry(i);
        // if (i % 1000 == 0) cout << i << "/" << nevt << endl;
        if(i % 10000 == 0) PrintMemoryUsage(Form("After processing %d events", i));

        //Redo the cluster to update the pi0 mass 
        for(int iblk = 0; iblk < 1080; iblk++) hasBlock[iblk] = false;

        Int_t nCaloBlock = caloev->GetNbBlocks();
        for(int iblk = 0; iblk < nCaloBlock; iblk++){ 
            TCaloBlock *block = caloev->GetBlock(iblk);
            Int_t nb = block->GetBlockNumber(); // Here the block number is in simulation numbering scheme
            Float_t energy = block->GetEnergy(0)*corr_pi0[nb];

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
        m_pi0_recons = (photon1 + photon2).M();

        if(photon1.E() > 1.4 && photon2.E() > 1.4) h_pi0M->Fill(m_pi0_recons);

        ev->Reset();
    }// event loop

    Double_t meanfit = h_pi0M->GetBinCenter(h_pi0M->GetMaximumBin());
    Double_t minfit = meanfit - 0.01;
    Double_t maxfit = meanfit + 0.01;
    TF1 *f_fit = new TF1(Form("f_fit_%d", iIter-1), "gausn(0)", minfit, maxfit);
    // f_fit->SetParLimits(0, 100, 100000);
    f_fit->SetParLimits(1, 0.12, 0.15);
    f_fit->SetParLimits(2, 0.00001, 0.03);
    // f_fit->SetParameters(10000, 0.13, 0.001, 1.2, 1.0); // N, mean, sigma, alpha, n
    // f_fit->SetParLimits(2, 0.00001, 0.03);
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
    for (Int_t i = 0; i < 1080; i++){
        log_Stream<<corr_pi0[i]<<endl;;
    }
    // log_Stream<<"pi0 calibration completed"<<endl;
    log_Stream<<"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
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