#include "/home/hhuang/workspace/group_h/analysis/MyHeader/Analysis.h"
#include "/group/nps/hhuang/software/NPS_SOFT/TDVCSDB.h"

void calibTree(int run_number, int iseg)
{   
    if((run_number == 4301 && iseg == 2) 
    || (run_number == 4301 && iseg == 4)
    || (run_number == 4302 && iseg == 0)
    || (run_number == 4303 && iseg == 3)
    || (run_number == 4359 && iseg == 0)
    || (run_number == 4361 && iseg == 4)
    || (run_number == 4399 && iseg == 0)
    || (run_number == 4408 && iseg == 1)){
        cout<<"Warning: can't find rootfile for run "<<run_number<<", segment "<<iseg<<" previously. Skip this job."<<endl;
        return;
    }
    // Choose if you want to do with hcana or waveform-fit data
    bool use_wf = false; // true for waveform-fit, false for hcana data

    if(!use_wf) cout<<"You are using hcana data for calibration"<<endl;
    else cout<<"You are using waveform-fit data for calibration"<<endl;
    
    // Check if the segments exist
    if(system(Form("test -f /group/nps/hhuang/analysis/DVCS_NPS2023/DVCS_analysis/pi0Calib_wf/WF_Kin_36_2/nps_production_%d_%d_4_wf.root", run_number, iseg))){ // 1 if file does not exist
        cout<<Form("/group/nps/hhuang/analysis/DVCS_NPS2023/DVCS_analysis/pi0Calib_wf/WF_Kin_36_2/nps_production_%d_%d_4_wf.root", run_number, iseg)<<" does not exist!"<<endl;
        return;
    }

    TFile *infile = TFile::Open(Form("/group/nps/hhuang/analysis/DVCS_NPS2023/DVCS_analysis/pi0Calib_wf/WF_Kin_36_2/nps_production_%d_%d_4_wf.root", run_number, iseg));
    
    // T tree__________________________________________________________
    TTree *t_T = (TTree*)infile->Get("T");
    if (!t_T){
        std::cerr << "No T tree found in the file.\n";
        // infile->Close();
        return;
    }
    
    Double_t runNb_T;
    Double_t evtNb_T;
    Double_t H_gtr_dp;
    Double_t H_gtr_ph;
    Double_t H_gtr_th;
    Double_t H_gtr_px;
    Double_t H_gtr_py;
    Double_t H_gtr_pz;
    Double_t H_react_x;
    Double_t H_react_y;
    Double_t H_react_z;
    Double_t H_hod_beta;
    Double_t H_cal_etottracknorm;
    Double_t H_cer_npeSum;

    Int_t Ndata_NPS_cal_fly_adcCounter;
    Double_t NPS_cal_fly_adcCounter[1080];
    Double_t NPS_cal_fly_adcSampPulseAmp[1080];
    Double_t NPS_cal_fly_adcSampPulseInt[1080];
    Double_t NPS_cal_fly_adcSampPulseTime[1080];
    Double_t NPS_cal_fly_adcSampPed[1080];

    Int_t Ndata_NPS_cal_vtpClusE;
    Double_t NPS_cal_vtpClusE[1080];
    Double_t NPS_cal_vtpClusTime[1080];
    Double_t NPS_cal_vtpClusX[1080];
    Double_t NPS_cal_vtpClusY[1080];

    //Event Level variables
    t_T->SetBranchAddress("g.runnum", &runNb_T);
    t_T->SetBranchAddress("g.evnum", &evtNb_T); // Global event number of T tree
    // HMS information
    t_T->SetBranchAddress("H.gtr.dp", &H_gtr_dp);
    t_T->SetBranchAddress("H.gtr.ph", &H_gtr_ph);
    t_T->SetBranchAddress("H.gtr.th", &H_gtr_th);
    t_T->SetBranchAddress("H.gtr.px", &H_gtr_px);
    t_T->SetBranchAddress("H.gtr.py", &H_gtr_py);
    t_T->SetBranchAddress("H.gtr.pz", &H_gtr_pz);
    t_T->SetBranchAddress("H.react.x", &H_react_x);
    t_T->SetBranchAddress("H.react.y", &H_react_y);
    t_T->SetBranchAddress("H.react.z", &H_react_z);
    t_T->SetBranchAddress("H.hod.beta", &H_hod_beta);
    t_T->SetBranchAddress("H.cal.etottracknorm", &H_cal_etottracknorm);
    t_T->SetBranchAddress("H.cer.npeSum", &H_cer_npeSum);

    // NPS fADC information
    t_T->SetBranchAddress("Ndata.NPS.cal.fly.adcCounter", &Ndata_NPS_cal_fly_adcCounter);
    t_T->SetBranchAddress("NPS.cal.fly.adcCounter", &NPS_cal_fly_adcCounter);
    t_T->SetBranchAddress("NPS.cal.fly.adcSampPulseAmp", &NPS_cal_fly_adcSampPulseAmp);
    t_T->SetBranchAddress("NPS.cal.fly.adcSampPulseInt", &NPS_cal_fly_adcSampPulseInt);
    t_T->SetBranchAddress("NPS.cal.fly.adcSampPulseTime", &NPS_cal_fly_adcSampPulseTime);
    t_T->SetBranchAddress("NPS.cal.fly.adcSampPed", &NPS_cal_fly_adcSampPed);

    // NPS VTP information
    t_T->SetBranchAddress("Ndata.NPS.cal.vtpClusE", &Ndata_NPS_cal_vtpClusE);
    t_T->SetBranchAddress("NPS.cal.vtpClusE", &NPS_cal_vtpClusE);
    t_T->SetBranchAddress("NPS.cal.vtpClusTime", &NPS_cal_vtpClusTime);
    t_T->SetBranchAddress("NPS.cal.vtpClusX", &NPS_cal_vtpClusX);
    t_T->SetBranchAddress("NPS.cal.vtpClusY", &NPS_cal_vtpClusY);

    // Disabla all and turn on the branches we need
    t_T->SetBranchStatus("*", false);
    t_T->SetBranchStatus("g.runnum", true);
    t_T->SetBranchStatus("g.evnum", true);
    t_T->SetBranchStatus("H.gtr.dp", true);
    t_T->SetBranchStatus("H.gtr.ph", true);
    t_T->SetBranchStatus("H.gtr.px", true);
    t_T->SetBranchStatus("H.gtr.py", true);
    t_T->SetBranchStatus("H.gtr.pz", true);
    t_T->SetBranchStatus("H.react.x", true);
    t_T->SetBranchStatus("H.react.y", true);
    t_T->SetBranchStatus("H.react.z", true);
    t_T->SetBranchStatus("H.hod.beta", true);
    t_T->SetBranchStatus("H.cal.etottracknorm", true);
    t_T->SetBranchStatus("H.cer.npeSum", true);

    // NPS fADC information
    t_T->SetBranchStatus("Ndata.NPS.cal.fly.adcCounter", true);
    t_T->SetBranchStatus("NPS.cal.fly.adcCounter", true);
    t_T->SetBranchStatus("NPS.cal.fly.adcSampPulseAmp", true);
    t_T->SetBranchStatus("NPS.cal.fly.adcSampPulseInt", true);
    t_T->SetBranchStatus("NPS.cal.fly.adcSampPulseTime", true);
    t_T->SetBranchStatus("NPS.cal.fly.adcSampPed", true);

    // NPS VTP information
    t_T->SetBranchStatus("Ndata.NPS.cal.vtpClusE", true);
    t_T->SetBranchStatus("NPS.cal.vtpClusE", true);
    t_T->SetBranchStatus("NPS.cal.vtpClusTime", true);
    t_T->SetBranchStatus("NPS.cal.vtpClusX", true);
    t_T->SetBranchStatus("NPS.cal.vtpClusY", true);

    // TChain of waveform tree___________________________________________________________
    TTree *t_wf = (TTree*)infile->Get("WF");
    if (!t_wf) {
        std::cerr << "No WF tree found in the file.\n";
        // infile->Close();
        return;
    }

    Double_t evtNb_wf;
    vector<Double_t> *chi2 = nullptr;
    vector<Int_t> *wfnpulse = nullptr;
    vector<Double_t> *wfampl = nullptr;
    vector<Double_t> *wftime = nullptr;
    t_wf->SetBranchAddress("evt", &evtNb_wf); // Global event number of wf tree
    t_wf->SetBranchAddress("chi2", &chi2); // Number of pulses found by TSpectrum in this block
    t_wf->SetBranchAddress("wfnpulse", &wfnpulse); // Number of pulses found by TSpectrum in this block 
    t_wf->SetBranchAddress("wfampl", &wfampl); // Waveform amplitude of each fitted pulse (mV)
    t_wf->SetBranchAddress("wftime", &wftime); // Pulse time of each fitted pulse per block (ns)

    // Use the TVirtualIndex to access the entry corresponding to the on in T tree
    TVirtualIndex *vIdx = t_wf->GetTreeIndex();
    TTreeIndex *idx = dynamic_cast<TTreeIndex *>(vIdx);
    if (!idx)
    {
        std::cerr << "No TTreeIndex found in the output tree.\n";
        infile->Close();
        return;
    }

    // Get the index array
    Long64_t *indexArray = idx->GetIndex();
    Int_t nevt_wf = t_wf->GetEntries();
    Int_t nevt_T = t_T->GetEntries();
    Double_t lastevent = 0.;

    if (nevt_wf != nevt_T){ // Check the number of events in both trees
        std::cerr << "Number of events in WF tree (" << nevt_wf << ") does not match number of events in T tree (" << nevt_T << ").\n";
        return;
    }
    
    // Timing offsets from Mark
    // ifstream fTimingOffset("/group/nps/mathison/analysis/NPS_Offsets/wfOffsets_3020.txt"); // x36_2
    // ifstream fTimingOffset("/group/nps/mathison/analysis/NPS_Offsets/wfOffsets_4253.txt"); // x60_4b
    // if(!fTimingOffset){
    //     cout<<"ERROR: Can't find the timing offsets file!!!"<<endl;
    //     return;
    // }
    // Double_t timingOffset[1080];
    // for(int iblk = 0; iblk < 1080; iblk++) fTimingOffset >> timingOffset[iblk];

    // For testing timing offsets==================

    // Timing offsets from Mark
    ifstream fTimingOffset("/group/nps/mathison/analysis/NPS_Offsets/wfOffsets_3020.txt"); // x36_2
    if(!fTimingOffset){
        cout<<"ERROR: Can't find the timing offsets file!!!"<<endl;
        return;
    }
    Double_t timingOffset_Mark[1080];
    for(int iblk = 0; iblk < 1080; iblk++) fTimingOffset >> timingOffset_Mark[iblk];

    // My timing offsets
    TFile *wfOffset_hao = TFile::Open("/group/nps/hhuang/analysis/DVCS_NPS2023/DVCS_analysis/pi0Calib_wf/checkCalib/Timing/plots/0p3GeV_v1/FitBlockTiming.root");
    TH1F *h_offset = (TH1F*)wfOffset_hao->Get(Form("h_fitTiming_noTimeCorr_run%d", run_number));
    Double_t timingOffset[1080];
    for(int iblk = 0; iblk < 1080; iblk++){
        if(!use_wf) timingOffset[iblk] = 0.;
        if(use_wf){
            Int_t icol = iblk%30;
            Int_t irow = (iblk-icol)/30;
            if((run_number == 3013 && (icol == 4 || iblk == 26))
            || (run_number == 3014 && (icol == 4)) 
            || (run_number == 3015 && (icol == 4 || iblk == 6 || iblk == 12 || iblk == 15 || iblk == 23 || iblk == 26 || iblk == 1056 || iblk == 1057 || iblk == 1089 || iblk == 1068 || iblk == 1070))
            || (run_number == 3016 && (icol == 4 || irow == 0 || irow == 35 || iblk == 65 || iblk == 126 || iblk == 245 || iblk == 1025 || iblk == 1026 || iblk == 1028))
            || (run_number == 3020 && (icol == 4 || irow == 35 || (0 <= iblk && iblk <= 26) || iblk == 52 || iblk == 1032 || (1050 <= iblk && iblk <= 1074)))) timingOffset[iblk] = timingOffset_Mark[iblk];
            else timingOffset[iblk] = -1*h_offset->GetBinContent(iblk+1);
        }
    }

    // Connect to the database and get kinematics variables________________________________________________
    gSystem->Load("/group/nps/hhuang/software/NPS_SOFT/libDVCS.so");
    TDVCSDB *db = new TDVCSDB("dvcs", "clrlpc", 3306, "hhuang", "");

    Double_t Beam_energy = *db->GetEntry_d("BEAM_param_Energy", run_number); // Beam energy in GeV
    Double_t HMS_mom = *db->GetEntry_d("SIMU_param_HMSmomentum", run_number); // HMS central momentum in GeV/c
    Double_t HMS_angle = *db->GetEntry_d("SIMU_param_HMSangle", run_number); // HMS angle in rad
    Double_t Target_amu = *db->GetEntry_d("TARGET_param_Amu", run_number); // Target amu

    // The angle and distance of calorimeter
    Double_t NPS_dist = *db->GetEntry_d("CALO_geom_Dist", run_number); // NPS distence in cm
    Double_t NPS_angle = *db->GetEntry_d("CALO_geom_Yaw", run_number); // NPS angle in rad
    Int_t *caloMaskBlock = new Int_t[1080]; // Get the mask block information in NPS numbering Scheme
    caloMaskBlock = db->GetEntry_i("CALO_flag_MaskBlock", run_number);

    // Elastic coefficients
    Double_t *coefElas = new Double_t[1080]; // Get the elastic coefficients (GeV/mV) in NPS numbering Scheme
    coefElas = db->GetEntry_d("CALO_calib_ElasCoef", run_number);
    for(int iblk = 0; iblk < 1080; iblk++){
        if(coefElas[iblk] <= 0 && !caloMaskBlock[iblk]) coefElas[iblk] = 0.014; // if the coefficient <= 0, set it to 0.014 GeV/mV
    }

    cout<<"Start clustering for Run "<<run_number<<", segment "<<iseg<<endl;
    cout<<"========== Run information =========="<<endl;
    cout<<"Beam energy: "<<Beam_energy<<" GeV"<<endl;
    cout<<"HMS momentum: "<<HMS_mom<<" GeV/c"<<endl;
    cout<<"HMS angle: "<<HMS_angle*TMath::RadToDeg()<<" degrees"<<endl;
    cout<<"Target amu: "<<Target_amu<<" GeV/c^2"<<endl;
    cout<<"NPS distance: "<<NPS_dist<<" cm"<<endl;
    cout<<"NPS angle: "<<NPS_angle*TMath::RadToDeg()<<" degrees"<<endl;
    cout<<"======================================"<<endl;
    cout<<endl;
    cout<<"========== NPS maske blocks =========="<<endl;
    for(int i = 0; i < 1080; i++) if(caloMaskBlock[i] == 1) cout<<i<<" ";
    cout<<endl;
    cout<<"======================================"<<endl;
    cout<<endl;
    cout<<"========== NPS elastic coefficients =========="<<endl;
    for(int i = 0; i < 1080; i++) cout<<coefElas[i]<<" ";
    cout<<endl;
    cout<<"======================================"<<endl;
    cout<<endl;
    cout<<"========== NPS timing offsets =========="<<endl;
    if(!use_wf) cout<<"Using hcana data, timing offsets are applied in the nps replay";
    if(use_wf) for(int iblk = 0; iblk < 1080; iblk++) cout<<timingOffset[iblk]<<" ";
    cout<<endl;
    cout<<"======================================"<<endl;
    cout<<endl;

    // Initial settings for clustering__________________________________________________
    TLorentzVector beam(0, 0, Beam_energy, Beam_energy);
    TLorentzVector p0(0, 0, 0, m_p);

    TDVCSEvent *ev = new TDVCSEvent();
    TCaloEvent *caloev = new TCaloEvent(run_number);
    ev->GetGeometry()->SetCaloTheta(-1*NPS_angle); // rad Note: the angle have to be negative here
    ev->GetGeometry()->SetCaloDist(NPS_dist); // cm

    // Create histograms__________________________________________________________
    TH1F *h_evtStat = new TH1F("h_evtStat", "Number of exclusive events;Run number;Counts", 3500, 1500.5, 5000.5);

    TH1F *h_reactX = new TH1F("h_reactX", "H.react.x;Reaction X [cm];Counts", 400, -0.4, 0.1);
    TH1F *h_reactY = new TH1F("h_reactY", "H.react.y;Reaction Y [cm];Counts", 400, -0.2, 0.2);
    TH1F *h_reactZ = new TH1F("h_reactZ", "H.react.z;Reaction Z [cm];Counts", 400, -20, 20);

    TH1F *h_th = new TH1F("h_th", "th", 100, -1, 1);
    TH2F *hh_dp_ph = new TH2F("hh_dp_ph", "dp vs. ph;", 100, -0.05, 0.05, 300, -15, 15);

    TH1F *h_etottracknorm = new TH1F("h_etottracknorm", "H.cal.etottracknorm", 1000, 0, 10);
    TH1F *h_npeSum = new TH1F("h_npeSum", "H.cer.npeSum", 1000, 0, 10);
    TH1F *h_beta = new TH1F("h_beta", "H.hod.beta;#beta;Counts", 1500, 0, 1.5);
    TH1F *h_NpsTime;
    if(!use_wf) h_NpsTime = new TH1F("h_NpsTime", ";Pulse time [ns];Number of events", 20000, 50, 250); // When using the T tree
    if(use_wf) h_NpsTime = new TH1F("h_NpsTime", "NPS timing with offset correction (ampl. > 10 mV);Pulse time [ns];Number of events", 20000, -100, 100); // When using the waveform tree

    TH1F *h_pi0M = new TH1F("h_pi0M", "M_{#gamma#gamma};M_{#gamma#gamma} [GeV];Counts", 150, 0.05, 0.2);
    h_pi0M->Sumw2();
    TH1F *h_Mx2 = new TH1F("h_Mx2", "#pi^{0} missing mass squire; M_{x}^{2} [GeV^{2}]; Counts", 180, 0, 1.8);
    h_Mx2->Sumw2();
    TH2F *hh_pi0M_Mx2 = new TH2F("hh_pi0M_Mx2", "M_{#gamma#gamma} vs. #pi^{0} missing mass squire;M_{#gamma#gamma} [GeV];M_{x}^{2} [GeV^{2}]", 200, 0.12, 0.16, 200, 0, 1.8);

    // Ontput files for Tree___________________________________________________________
    TFile *outfile;
    if(!use_wf) outfile = new TFile(Form("/group/nps/hhuang/analysis/DVCS_NPS2023/DVCS_analysis/pi0Calib_wf/calibTree/x36_2_1/prodTree_pass2_v3_%d_%d.root", run_number, iseg), "recreate"); // When using the T tree
    if(use_wf) outfile = new TFile(Form("/group/nps/hhuang/analysis/DVCS_NPS2023/DVCS_analysis/pi0Calib_wf/calibTree/x36_2_1/prodTree_wf_v7_%d_%d.root", run_number, iseg), "recreate");// When using the waveform tree

    if (!outfile){ // Check if the output file was created successfully
        std::cerr << "Error creating output file.\n";
        return;
    }

    TTree *t_prod = new TTree("t_prod", "t_prod");
    // HMS information
    Double_t hms_Vx;
    Double_t hms_Vy;
    Double_t hms_Vz;
    t_prod->Branch("hms_Vx", &hms_Vx);
    t_prod->Branch("hms_Vy", &hms_Vy);
    t_prod->Branch("hms_Vz", &hms_Vz);
    // NPS information
    // Double_t clust_ene1;
    // Double_t clust_ene2;
    t_prod->Branch("caloev", &caloev);
    // t_prod->Branch("clust_ene1", &clust_ene1);
    // t_prod->Branch("clust_ene2", &clust_ene2);

    // Event loop
    Int_t nevt = nevt_T;
    // Int_t nevt = 100000; // for testing

    Int_t nExcEvt = 0;
    for(Int_t ievt = 0; ievt < nevt; ievt++){
        if(run_number == 4365 && iseg == 3 && ievt == 233096) continue; // There is some problem with this event in x60_4b, skip it
        t_T->GetEntry(ievt);
        t_wf->GetEntry(indexArray[ievt]); // indexArray[i] is the original-entry number for the i-th smallest evt
        if(ievt%100000==0) cout << "Looking at entry = " << ievt << "  (" << 100.*ievt/nevt_T << "%)" << endl;

        // Fill the histograms before selection
        h_reactX->Fill(H_react_x);
        h_reactY->Fill(H_react_y);
        h_reactZ->Fill(H_react_z);
        h_th->Fill(H_gtr_th);
        hh_dp_ph->Fill(H_gtr_ph, H_gtr_dp);
        h_etottracknorm->Fill(H_cal_etottracknorm);
        h_npeSum->Fill(H_cer_npeSum);
        h_beta->Fill(H_hod_beta);

        // Event selection for pi0 reconstruction
        if(!(TMath::Abs(H_react_z) < 4)) continue;
        if(!(TMath::Abs(H_gtr_ph) < 0.03 && TMath::Abs(H_gtr_dp) < 10)) continue;
        if(!(TMath::Abs(H_hod_beta-1.) < 0.2)) continue;
        if(!(H_cer_npeSum > 2)) continue;
        if(!(H_cal_etottracknorm > 0.6)) continue;

        
        if(!use_wf){ // When using the T tree____________________________________________________
            Bool_t hasPulse[nblk];
            for(int iblk = 0; iblk < nblk; iblk++) hasPulse[iblk] = false;

            Int_t Ndata = Ndata_NPS_cal_fly_adcCounter;
            for(int idata = 0; idata < Ndata; idata++){
                if(!(0 <= NPS_cal_fly_adcCounter[idata])) continue;
                if(!(NPS_cal_fly_adcCounter[idata] < 1080)) continue;
                if(!(abs(NPS_cal_fly_adcSampPulseTime[idata] - 150) < 3)) continue;
                
                Int_t iblk = NPS_cal_fly_adcCounter[idata]; // block number with current numbering scheme
                Int_t iblk_sim = bnConv_NewToOld(iblk); // block number with simulation numbering scheme
                TCaloBlock *block = caloev->AddBlock(iblk_sim); // add block with simulation numbering scheme
                block->AddPulse(NPS_cal_fly_adcSampPulseAmp[idata] * coefElas[iblk], 0); // Energy,Time
                hasPulse[iblk_sim] = true; // simulation numbering scheme
            }

            for(int iblk_sim = 0; iblk_sim < nblk; iblk_sim++){
                if(!hasPulse[iblk_sim]) caloev->AddBlock(iblk_sim); // simulation numbering scheme
            }
        }
        

        if(use_wf){ // When using the waveform tree_____________________________________________
            Int_t flatIdx = 0; // accumulated number of pulses in previous blocks
            for (size_t iblk=0; iblk<wfnpulse->size(); iblk++){
                Int_t nPulse = (*wfnpulse)[iblk];
                Int_t iblk_sim = bnConv_NewToOld(iblk); // block number with simulation numbering scheme
                TCaloBlock *block = caloev->AddBlock(iblk_sim); // add block with simulation numbering scheme
                if (nPulse > 0 && (*chi2)[iblk] > 0 && !caloMaskBlock[iblk]){ // if we want to read the pulse in this block
                    Int_t pulseIdx = flatIdx + 0; // index of the first pulse in this block
                    for (Int_t p = 0; p < nPulse; p++){ // loop over the pulses in this block
                        if((*wfampl)[flatIdx + p] > 10) h_NpsTime->Fill((*wftime)[flatIdx + p]+timingOffset[iblk]); // Fill the histogram with timing offset correction
                    }
                    if(nPulse > 1){ // Look at the pulse closest to the coincidence time
                        for (Int_t p = 0; p < nPulse; p++){ // loop over the pulses in this block
                            if(abs((*wftime)[flatIdx + p]+timingOffset[iblk]) <= abs((*wftime)[pulseIdx]+timingOffset[iblk])) pulseIdx = flatIdx + p; // Find the pulse closest to the coincidence time
                        }
                    }
                    if(abs((*wftime)[pulseIdx]+timingOffset[iblk]) < 3) block->AddPulse((*wfampl)[pulseIdx] * coefElas[iblk], 0 /*(*wftime)[pulseIdx]+timingOffset[iblk]*/); // Add energy and timing to the block
                } // block selection
                flatIdx += nPulse;
            } // loop over 1080 blocks
        }

        caloev->TriggerSim(0.1);     // Energy threshold 0.05 GeV of every 2x2 for clustering
        caloev->DoClustering(-3, 3); // Time window (-3,3) [ns]
        ev->SetCaloEvent(caloev);
        ev->SetVertex(H_react_x, H_react_y, H_react_z);
        // ev->SetVertex(0., 0., vertex * 100.);
        // clust_ene1 = -1;
        // clust_ene2 = -1;
        if(caloev->GetNbClusters() == 2){
            caloev->GetCluster(0)->Analyze(1);
            caloev->GetCluster(1)->Analyze(1);

            TLorentzVector photon1 = ev->GetPhoton(0, 7, 0);
            TLorentzVector photon2 = ev->GetPhoton(1, 7, 0);
            TLorentzVector kpvec;
            kpvec.SetXYZM(H_gtr_px, H_gtr_py, H_gtr_pz, 0.000511);
            Double_t pi0_m = (photon1 + photon2).M();
            Double_t pi0_mm2 = (beam + p0 - kpvec - photon1 - photon2).Mag2();

            // if (!(0.08 < pi0_m && pi0_m < 0.2 && pi0_mm2 < 1.2 && photon1.E() > 1 && photon2.E() > 1)) continue;
            if (!(0.05 < pi0_m && pi0_m < 0.2/* && photon1.E() > 1 && photon2.E() > 1*/)) continue;
            
            nExcEvt += 1;

            h_pi0M->Fill(pi0_m);
            h_Mx2->Fill(pi0_mm2);
            hh_pi0M_Mx2->Fill(pi0_m, pi0_mm2);

            hms_Vx = H_react_x;
            hms_Vy = H_react_y;
            hms_Vz = H_react_z;
            // clust_ene1 = caloev->GetCluster(0)->GetEnergy();
            // clust_ene2 = caloev->GetCluster(1)->GetEnergy();
            // hms_Px = H_gtr_px;
            // hms_Py = H_gtr_py;
            // hms_Pz = H_gtr_pz;
            t_prod->Fill();
        } // clus==2

        // reinitialization for next event
        caloev->Reset();
    } // Event loop

    h_evtStat->SetBinContent(run_number-1500, nExcEvt);

    outfile->cd();
    t_prod->Write();
    
    h_evtStat->Write();
    h_pi0M->Write();
    h_Mx2->Write();
    hh_pi0M_Mx2->Write();

    h_reactX->Write();
    h_reactY->Write();
    h_reactZ->Write();
    h_beta->Write();
    h_etottracknorm->Write();
    h_npeSum->Write();
    h_th->Write();
    hh_dp_ph->Write();
    h_NpsTime->Write();
    
    outfile->Close();
}
