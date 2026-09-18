#include "calibHeader/MyDB.h"
#include "calibHeader/DrawSetting.h"
#include "calibHeader/Analysis.h"

double extractF(const std::string &filename) {
    std::ifstream infile(filename);
    if(!infile.is_open()) {
        std::cerr << "Error: cannot open file " << filename << std::endl;
        return NAN;
    }

    std::string line;
    while(std::getline(infile, line)) {
        if (line.find("F =") != std::string::npos) {
            std::stringstream ss(line);
            std::string eq, equal_sign;
            double F_value;
            ss >> eq >> equal_sign >> F_value;
            return F_value;
        }
    }

    std::cerr << "Warning: F not found in " << filename << std::endl;
    return NAN;
}

void Draw_pi0Calib(TString Kine, int TargetFlag, int startcycle, int group, bool doCyclePlot = false)
{
    Global_setting();
    TGaxis::SetMaxDigits(4);
    TH1::SetDefaultSumw2();

    const Int_t nIter = 8;
    const Int_t nTest = group; // number of cycles to check statistics
    
    TString Tar;
    if(TargetFlag == 0) Tar = "LH2";
    else if(TargetFlag == 1) Tar = "LD2";
    else if(TargetFlag == -1) Tar = "LH2_LD2";
    else{
        cout<<"ERROR: Unknown target!!!"<<endl;
        return;
    }

    // Folder name of calibration results
    TString prefix_1 = "wf";
    // TString prefix_1 = "pass2";
    TString prefix_2 = Form("cycle%d", startcycle);
    for(int ic = startcycle+1; ic < startcycle+group; ic++) prefix_2+=Form("_cycle%d", ic);
    TString filename = Form("x%s_%s_%s_%s", Kine.Data(), Tar.Data(), prefix_1.Data(), prefix_2.Data());

    // make folder for the plots
    system(Form("mkdir Result/%s/plots", filename.Data()));

    // Get the list of runs
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

    // Coefficient histograms from calibration
    TFile *infile_summary = TFile::Open(Form("Result/%s/CalibSummary.root", filename.Data()));
    
    TH1F *h_coEff = (TH1F*)infile_summary->Get("h_coEff");
    TH2F *hh_coEff = (TH2F*)infile_summary->Get("hh_coEff");

    TH2F *hh_iblk = new TH2F("hh_iblk", "", 30, -0.5, 29.5, 36, -0.5, 35.5);

    for(int iblk = 0; iblk < 1080; iblk++){
        Int_t icol = iblk%30;
        Int_t irow = (iblk-icol)/30;
        hh_iblk->SetBinContent(icol+1, irow+1, iblk);
    }

    // pi0 invariant mass, missing mass, and energy resolution from Tree
    TTree *t_checkCalib = (TTree*)infile_summary->Get("t_checkCalib");

    // Cut to apply
    Double_t E_rmCorr = 20.;
    TString cut_accpt = "abs(clust_x1)<29.025 && abs(clust_x2)<29.025 && abs(clust_y1)<35.475 && abs(clust_y2)<35.475";
    TString cut_photE = "clust_ene1>1. && clust_ene2>1.";
    TString cut_exclusive = Form("0.13<pi0_m && pi0_m<0.14 && 0.8<pi0_mm2+%f*(pi0_m-0.1349766) && pi0_mm2+%f*(pi0_m-0.1349766)<1", E_rmCorr, E_rmCorr);

    // Histograms from the tree
    TH1F *h_m = new TH1F("h_m", "M_{#gamma#gamma};M_{#gamma#gamma} [GeV];Counts", 100, 0.08, 0.2);
    TH1F *h_mx2 = new TH1F("h_mx2", "M_{x}^{2};M_{x}^{2} [GeV^2];Counts", 100, 0, 1.5);
    TH1F *h_mx2_rmCorr = new TH1F("h_mx2_rmCorr", Form("M_{x}^{2} (correlation removed);M_{x}^{2}+%.0f #times (M_{#gamma#gamma}-M_{#pi^{0}}) [GeV^2];Counts", E_rmCorr), 100, 0, 1.5);

    TH2F *hh_m_vs_mx2 = new TH2F("hh_m_vs_mx2", "M_{#gamma#gamma} vs. M_{x}^{2};M_{x}^{2} [GeV^2];M_{#gamma#gamma} [GeV]", 100, 0, 1.5, 100, 0.08, 0.2);
    TH2F *hh_m_vs_mx2_rmCorr = new TH2F("hh_m_vs_mx2_rmCorr", Form("M_{#gamma#gamma} vs. M_{x}^{2} (correlation removed);M_{x}^{2}+%.0f#times(M_{#gamma#gamma}-M_{#pi^{0}}) [GeV^2];M_{#gamma#gamma} [GeV]", E_rmCorr), 100, 0, 1.5, 100, 0.08, 0.2);

    TH1F *h_diffEpi0Ecal = new TH1F("h_diffEpi0Ecal", "#DeltaE (E_{#gamma#gamma}-E_{#pi^{0}, predicted}, exclusive #pi^{0});E_{#gamma#gamma}-E_{#pi^{0}, predicted} [GeV]", 50, -1, 1);
    TH1F *h_Epi0Reso = new TH1F("h_Epi0Reso", "Energy resolution of exclusive #pi^{0};(E_{#gamma#gamma}-E_{#pi^{0}, predicted})/E_{#pi^{0}, predicted}", 50, -0.1, 0.1);

    t_checkCalib->Draw("pi0_m>>h_m", Form("%s && %s", cut_accpt.Data(), cut_photE.Data()),"colz");
    t_checkCalib->Draw("pi0_mm2>>h_mx2", Form("%s && %s", cut_accpt.Data(), cut_photE.Data()),"colz");
    t_checkCalib->Draw(Form("pi0_mm2+%f*(pi0_m-0.1349766)>>h_mx2_rmCorr", E_rmCorr), Form("%s && %s", cut_accpt.Data(), cut_photE.Data()),"colz");
    
    t_checkCalib->Draw("pi0_m:pi0_mm2>>hh_m_vs_mx2", Form("%s && %s", cut_accpt.Data(), cut_photE.Data()),"colz");
    t_checkCalib->Draw(Form("pi0_m:pi0_mm2+%f*(pi0_m-0.1349766)>>hh_m_vs_mx2_rmCorr", E_rmCorr), Form("%s && %s", cut_accpt.Data(), cut_photE.Data()),"colz");
    
    t_checkCalib->Draw("(clust_ene1+clust_ene2)-pi0_Ecal>>h_diffEpi0Ecal", Form("%s && %s && %s", cut_accpt.Data(), cut_exclusive.Data(), cut_photE.Data()),"colz");
    t_checkCalib->Draw("((clust_ene1+clust_ene2)-pi0_Ecal)/pi0_Ecal>>h_Epi0Reso", Form("%s && %s && %s", cut_accpt.Data(), cut_exclusive.Data(), cut_photE.Data()),"colz");

    TCanvas *c_tree = new TCanvas("c_tree", "", 1000, 1000);
    c_tree->SetGridx();

    format_Line(h_m, 1, kBlue, 3);
    format_Mark(h_m, 20, kBlue, 1.5);
    h_m->Draw("lpe");
    TLine *l_pi0 = new TLine(0.1349766, 0, 0.1349766, h_m->GetMaximum());
    format_Line(l_pi0, 7, kRed, 3);
    l_pi0->Draw("l same");
    c_tree->Print(Form("Result/%s/plots/CalibQuality.pdf(", filename.Data()), "Title:pi0M_summary");
    c_tree->Clear();

    format_Line(h_mx2, 1, kBlue, 3);
    format_Mark(h_mx2, 20, kBlue, 1.5);
    h_mx2->Draw("lpe");
    TLine *l_pi0_mx2 = new TLine(0.88, 0, 0.88, h_mx2->GetMaximum());
    format_Line(l_pi0_mx2, 7, kRed, 3);
    l_pi0_mx2->Draw("l same");
    c_tree->Print(Form("Result/%s/plots/CalibQuality.pdf", filename.Data()), "Title:pi0Mx2_summary");
    c_tree->Clear();

    format_Line(h_mx2_rmCorr, 1, kBlue, 3);
    format_Mark(h_mx2_rmCorr, 20, kBlue, 1.5);
    h_mx2_rmCorr->Draw("lpe");
    l_pi0_mx2->SetY2(h_mx2_rmCorr->GetMaximum());
    l_pi0_mx2->Draw("l same");
    c_tree->Print(Form("Result/%s/plots/CalibQuality.pdf", filename.Data()), "Title:pi0Mx2_rmCorr_summary");
    c_tree->Clear();

    format_Line(h_diffEpi0Ecal, 1, kBlue, 3);
    format_Mark(h_diffEpi0Ecal, 20, kBlue, 1.5);
    h_diffEpi0Ecal->Draw("lpe");
    c_tree->Print(Form("Result/%s/plots/CalibQuality.pdf", filename.Data()), "Title:pi0Ediff_summary");
    c_tree->Clear();

    format_Line(h_Epi0Reso, 1, kBlue, 3);
    format_Mark(h_Epi0Reso, 20, kBlue, 1.5);
    h_Epi0Reso->Draw("lpe");
    c_tree->Print(Form("Result/%s/plots/CalibQuality.pdf", filename.Data()), "Title:pi0Ereso_summary");
    c_tree->Clear();

    c_tree->SetGridy();
    format_Line(hh_m_vs_mx2, 1, kBlue, 3);
    format_Mark(hh_m_vs_mx2, 20, kBlue, 1.5);
    hh_m_vs_mx2->Draw("colz");
    c_tree->Print(Form("Result/%s/plots/CalibQuality.pdf", filename.Data()), "Title:pi0M_Mx2_summary");
    c_tree->Clear();

    format_Line(hh_m_vs_mx2_rmCorr, 1, kBlue, 3);
    format_Mark(hh_m_vs_mx2_rmCorr, 20, kBlue, 1.5);
    hh_m_vs_mx2_rmCorr->Draw("colz");
    c_tree->Print(Form("Result/%s/plots/CalibQuality.pdf)", filename.Data()), "Title:pi0M_Mx2_rmCorr_summary");
    c_tree->Clear();


    // Summary plot as a function of iterations
    TFile *infile[nIter+1];
    TH1 *h_pi0M;
    TF1 *f_fit;
    TLegend *leg_pi0M;
    TCanvas *c_pi0M_duringCalib = new TCanvas("c_pi0M_duringCalib", "", 1000, 1000);
    TCanvas *c_pi0M = new TCanvas("c_pi0M", "", 1000, 1000);
    TH1F *h_mass_iter = new TH1F("h_mass_iter", "M_{#pi^{0}} vs. iteration;Iteration;M_{#pi^{0}} [MeV/c^{2}]",nIter+1, -0.5, nIter+0.5);
    TH1F *h_width_iter = new TH1F("h_width_iter", "#sigma_{#pi^{0}} vs. iteration;Iteration;#sigma_{#pi^{0}} [MeV/c^{2}]",nIter+1, -0.5, nIter+0.5);
    TH1F *h_dmdpi0m_iter = new TH1F("h_dmdpi0m_iter", "|#DeltaM_{#pi^{0}}|/#deltaM_{fit} vs. iteration;Iteration;|#deltaM_{#pi^{0}}|/#sigma_{#pi^{0}}", nIter+1, -0.5, nIter+0.5);
    
    for(int iIter = 0; iIter < nIter+1; iIter++){
        infile[iIter] = TFile::Open(Form("Result/%s/f_Mgg_%d.root", filename.Data(), iIter));
        h_pi0M = (TH1F*)infile[iIter]->Get(Form("h_pi0M_%d", iIter));
        c_pi0M_duringCalib->cd();
        h_pi0M->Draw("lpe");
        if(iIter == 0) c_pi0M_duringCalib->Print(Form("Result/%s/plots/pi0M_duringCalib.pdf(", filename.Data()), Form("Title:pi0M_iter_%d", iIter));
        else if(iIter == nIter) c_pi0M_duringCalib->Print(Form("Result/%s/plots/pi0M_duringCalib.pdf)", filename.Data()), Form("Title:pi0M_iter_%d", iIter));
        else c_pi0M_duringCalib->Print(Form("Result/%s/plots/pi0M_duringCalib.pdf", filename.Data()), Form("Title:pi0M_iter_%d", iIter));

        h_pi0M->Rebin(2);
        Double_t meanfit = h_pi0M->GetBinCenter(h_pi0M->GetMaximumBin());
        Double_t minfit = meanfit - 0.035;
        Double_t maxfit = meanfit + 0.035;
        f_fit = new TF1(Form("f_fit_%d", iIter), "gausn(0)+pol1(3)", minfit, maxfit);

        f_fit->SetParLimits(0, 1, 10000);
        f_fit->SetParLimits(1, 0.1, 0.15);
        f_fit->SetParLimits(2, 0.0005, 0.015);
        h_pi0M->Fit(f_fit, "R");
        Double_t wBin = h_pi0M->GetBinWidth(1);
        Double_t N = f_fit->GetParameter(0)/wBin;
        Double_t N_err = f_fit->GetParError(0)/wBin;
        Double_t mean = f_fit->GetParameter(1)*1000;
        Double_t mean_err = f_fit->GetParError(1)*1000;
        Double_t sigm = f_fit->GetParameter(2)*1000;
        Double_t sigm_err = f_fit->GetParError(2)*1000;
        Double_t chi2 = f_fit->GetChisquare();
        Double_t ndf = f_fit->GetNDF();

        h_mass_iter->SetBinContent(iIter+1, mean);
        h_mass_iter->SetBinError(iIter+1, mean_err);
        h_width_iter->SetBinContent(iIter+1, sigm);
        h_width_iter->SetBinError(iIter+1, sigm_err);
        h_dmdpi0m_iter->SetBinContent(iIter+1, abs(mean-m_pi0*1000)/mean_err);
        cout<<mean<<", "<<m_pi0*1000<<", "<<mean_err<<", "<<abs(mean-m_pi0*1000)/mean_err<<endl;
        h_dmdpi0m_iter->SetBinError(iIter+1, 0.0001);

        leg_pi0M = new TLegend(0.40, 0.65, 0.85, 0.90);

        c_pi0M->cd();
        h_pi0M->GetYaxis()->SetRangeUser(0, 1.8*h_pi0M->GetMaximum());
        format_Mark(h_pi0M, 20, kBlack, 1.5);
        format_Line(h_pi0M, 1, kBlack, 3);
        h_pi0M->Draw("lpe");

        // Draw the pi0 mass before calibration
        if(iIter == nIter){
            TH1F *h_pi0M_0 = (TH1F*)infile[0]->Get("h_pi0M_0");
            // h_pi0M_0->Rebin(2);
            leg_pi0M->AddEntry(h_pi0M_0, "Before calibration", "l");
            leg_pi0M->AddEntry(h_pi0M, "After calibration", "lpe");
            leg_pi0M->AddEntry(f_fit, "Fit (gausn+pol1)", "l");
            leg_pi0M->AddEntry((TObject*)0, Form("#chi^{2}/NDF = %.2f", chi2/ndf), "");
            leg_pi0M->AddEntry((TObject*)0, Form("N_{#pi^{0}} = %.0f #pm %.0f", N, N_err), "");
            leg_pi0M->AddEntry((TObject*)0, Form("M_{#pi^{0}} = %.2f #pm %.2f MeV", mean, mean_err), "");
            leg_pi0M->AddEntry((TObject*)0, Form("#sigma_{#pi^{0}} = %.2f #pm %.2f MeV", sigm, sigm_err), "");
            leg_pi0M->SetY1(0.6);
            format_Line(h_pi0M_0, 1, kGreen+2, 3);
            h_pi0M_0->Draw("same hist");
        }

        else{
            leg_pi0M->AddEntry(h_pi0M, "Data", "lpe");
            leg_pi0M->AddEntry(f_fit, "Fit (gausn+pol1)", "l");
            leg_pi0M->AddEntry((TObject*)0, Form("#chi^{2}/NDF = %.2f", chi2/ndf), "");
            leg_pi0M->AddEntry((TObject*)0, Form("N_{#pi^{0}} = %.0f #pm %.0f", N, N_err), "");
            leg_pi0M->AddEntry((TObject*)0, Form("M_{#pi^{0}} = %.2f #pm %.2f MeV", mean, mean_err), "");
            leg_pi0M->AddEntry((TObject*)0, Form("#sigma_{#pi^{0}} = %.2f #pm %.2f MeV", sigm, sigm_err), "");
        }

        format_Line(f_fit, 1, kRed, 5);
        f_fit->Draw("C same");

        format_Legend(leg_pi0M);
        leg_pi0M->Draw();

        if(iIter == 0) c_pi0M->Print(Form("Result/%s/plots/pi0M_afterCalib.pdf(", filename.Data()), Form("Title:pi0M_iter_%d", iIter));
        else c_pi0M->Print(Form("Result/%s/plots/pi0M_afterCalib.pdf", filename.Data()), Form("Title:pi0M_iter_%d", iIter));
        c_pi0M->Clear();
    }

    TLine *l_mass = new TLine(-0.5, 1000*0.1349766, nIter+0.5, 1000*0.1349766);

    TLegend *leg_mass = new TLegend(0.45, 0.80, 0.85, 0.85);
    // leg_mass->SetHeader("#bf{E_{photons} > 1.3 GeV}");
    leg_mass->AddEntry(l_mass, "M_{#pi^{0}} = 0.1349766 GeV/c^{2}", "l");

    // TLegend *leg_width = new TLegend(0.45, 0.80, 0.85, 0.85);
    // leg_width->SetHeader("#bf{E_{photons} > 1.3 GeV}");

    c_pi0M->cd();
    h_mass_iter->GetYaxis()->SetRangeUser(h_mass_iter->GetMinimum()-5, h_mass_iter->GetMaximum()+5);
    format_Mark(h_mass_iter, 20, kBlue, 1.5);
    format_Line(h_mass_iter, 1, kBlue, 3);
    h_mass_iter->Draw("lpe");

    format_Line(l_mass, 7, kRed, 3);
    l_mass->Draw("l same");

    format_Legend(leg_mass);
    leg_mass->Draw();

    c_pi0M->Print(Form("Result/%s/plots/pi0M_afterCalib.pdf", filename.Data()), "Title:Mass_iter");
    c_pi0M->Clear();

    h_width_iter->GetYaxis()->SetRangeUser(0.8*h_width_iter->GetMinimum(), 1.2*h_width_iter->GetMaximum());
    format_Mark(h_width_iter, 20, kBlue, 1.5);
    format_Line(h_width_iter, 1, kBlue, 3);
    h_width_iter->Draw("lpe");

    // format_Legend(leg_width);
    // leg_width->Draw();
    c_pi0M->Print(Form("Result/%s/plots/pi0M_afterCalib.pdf", filename.Data()), "Title:Width_iter");

    // Additional correction after calibration
    TH1F *h_pi0MeanvsRun_old = (TH1F*)infile_summary->Get("h_pi0MeanvsRun_old");
    TH1F *h_pi0MeanvsRun_new = (TH1F*)infile_summary->Get("h_pi0MeanvsRun_new");

    TH1F *h_pi0SigmvsRun_old = (TH1F*)infile_summary->Get("h_pi0SigmvsRun_old");
    TH1F *h_pi0SigmvsRun_new = (TH1F*)infile_summary->Get("h_pi0SigmvsRun_new");

    h_pi0MeanvsRun_old->SetTitle(Form("M_{#pi^{0}} after additional correction;Run number;M_{#pi^{0}} [MeV/c^{2}]"));
    h_pi0MeanvsRun_old->GetXaxis()->SetLabelSize(0.03);
    h_pi0SigmvsRun_old->GetXaxis()->SetLabelSize(0.03);

    h_pi0SigmvsRun_old->SetTitle(Form("#sigma_{#pi^{0}} after additional correction;Run number;#sigma_{#pi^{0}} [MeV/c^{2}]"));

    TLegend *leg_addCorr = new TLegend(0.15, 0.15, 0.45, 0.25);
    leg_addCorr->AddEntry(h_pi0MeanvsRun_old, "Before additional correction", "lpe");
    leg_addCorr->AddEntry(h_pi0MeanvsRun_new, "After additional correction", "lpe");

    c_pi0M->cd();
    c_pi0M->Clear();
    c_pi0M->SetGridy();

    h_pi0MeanvsRun_old->GetYaxis()->SetRangeUser(100, 150);
    format_Mark(h_pi0MeanvsRun_old, 24, kRed, 1.5);
    format_Line(h_pi0MeanvsRun_old, 1, kRed, 3);
    h_pi0MeanvsRun_old->Draw("lpe");
    format_Mark(h_pi0MeanvsRun_new, 24, kBlue, 1.5);
    format_Line(h_pi0MeanvsRun_new, 1, kBlue, 3);
    h_pi0MeanvsRun_new->Draw("lpe same");

    format_Legend(leg_addCorr);
    leg_addCorr->Draw();
    
    c_pi0M->Print(Form("Result/%s/plots/pi0M_afterCalib.pdf", filename.Data()), "Title:Mean_additionalCorr");

    c_pi0M->cd();
    c_pi0M->Clear();
    c_pi0M->SetGridy();

    h_pi0SigmvsRun_old->GetYaxis()->SetRangeUser(2, 6);
    format_Mark(h_pi0SigmvsRun_old, 24, kRed, 1.5);
    format_Line(h_pi0SigmvsRun_old, 1, kRed, 3);
    h_pi0SigmvsRun_old->Draw("lpe");
    format_Mark(h_pi0SigmvsRun_new, 24, kBlue, 1.5);
    format_Line(h_pi0SigmvsRun_new, 1, kBlue, 3);
    h_pi0SigmvsRun_new->Draw("lpe same");

    leg_addCorr->Draw();

    c_pi0M->Print(Form("Result/%s/plots/pi0M_afterCalib.pdf)", filename.Data()), "Title:Width_additionalCorr");

    TCanvas *c_coEff = new TCanvas("c_coEff", "", 1000, 1000);
    format_Line(h_coEff, 1, kBlue, 3);
    h_coEff->Draw("lpe");
    c_coEff->Print(Form("Result/%s/plots/SummyAndDebugging.pdf(", filename.Data()), "Title:coefficient");
    
    h_coEff->GetYaxis()->SetRangeUser(0, 2);
    // c_coEff->Print(Form("Result/%s/plots/SummyAndDebugging.pdf", filename.Data()));

    TCanvas *c_coEff2D = new TCanvas("c_coEff2D", "", 2000, 1000);
    c_coEff2D->SetLeftMargin(5);
    hh_coEff->Draw("colz");
    hh_iblk->Draw("text same");
    c_coEff2D->Print(Form("Result/%s/plots/SummyAndDebugging.pdf", filename.Data()), "Title:coefficient2D");

    hh_coEff->GetZaxis()->SetRangeUser(0, 2);
    hh_iblk->Draw("text same");
    // c_coEff2D->Print(Form("Result/%s/plots/SummyAndDebugging.pdf", filename.Data()));

    TCanvas *c_dm = new TCanvas("c_dm", "", 1000, 1000);
    h_dmdpi0m_iter->GetXaxis()->SetRangeUser(0.5, nIter+0.5);
    h_dmdpi0m_iter->GetYaxis()->SetRangeUser(0, 15);
    format_Mark(h_dmdpi0m_iter, 20, kBlue, 1.5);
    format_Line(h_dmdpi0m_iter, 1, kBlue, 3);
    h_dmdpi0m_iter->Draw("lpe");
    c_dm->Print(Form("Result/%s/plots/SummyAndDebugging.pdf", filename.Data()), "Title:dmdsigma_iter");

    // epsilon and F for debugging
    TFile *infile_debug[nIter];
    TH1F *h_epsilon_simBn[nIter];
    Double_t epsilon_simBn[1080];
    TH1F *h_epsilon[nIter];
    TH2F *hh_epsilon[nIter];
    Double_t epsilon_avg[nIter];
    Double_t F_value[nIter];
    TH2D *hh_mat[nIter], *hh_mat_inv[nIter], *hh_mat_mult[nIter];
    TH1F *h_epsilon_iter = new TH1F("h_epsilon_iter", "Average #varepsilon vs. iteration;Iteration;<#varepsilon>", nIter, 0.5, nIter+0.5);
    // TH1F *h_F_iter = new TH1F("h_F_iter", "F vs. iteration;Iteration;F", nIter, 0.5, nIter+0.5);
    for(int i = 0; i<nIter; i++){
        epsilon_avg[i] = 0;
        for(int iblk = 0; iblk<1080; iblk++) epsilon_simBn[iblk] = 0;
        F_value[i] = 0;

        infile_debug[i] = TFile::Open(Form("Result/%s/f_Debug_%d.root", filename.Data(), i+1));
        h_epsilon_simBn[i] = (TH1F*)infile_debug[i]->Get(Form("h_epsilon_%d", i+1));
        h_epsilon[i] = new TH1F(Form("h_epsilon_nps_%d", i+1), "#varepsilon of #pi^{0} calibration;Block number;#varepsilon", 1080, -0.5, 1079.5);
        hh_epsilon[i] = new TH2F(Form("hh_epsilon_%d", i+1), "#varepsilon of #pi^{0} calibration;Column number, Row number", 30, -0.5, 29.5, 36, -0.5, 35.5);
        for(int iblk = 0; iblk<1080; iblk++) epsilon_simBn[iblk] = h_epsilon_simBn[i]->GetBinContent(iblk+1); // epsilon in simulation numbering scheme
        for(int iblk = 0; iblk<1080; iblk++){ // NPS numbering scheme
            Int_t iblk_sim = bnConv_NewToOld(iblk);
            Int_t icol = iblk%30;
            Int_t irow = (iblk-icol)/30;
            h_epsilon[i]->SetBinContent(iblk+1, epsilon_simBn[iblk_sim]);
            hh_epsilon[i]->SetBinContent(icol+1, irow+1, epsilon_simBn[iblk_sim]);
            epsilon_avg[i] += epsilon_simBn[iblk_sim];
        }
        epsilon_avg[i] /= 1080;
        h_epsilon_iter->SetBinContent(i+1, epsilon_avg[i]);
        h_epsilon_iter->SetBinError(i+1, 0.00000000001);

        // F value
        // std::string filename = "Result/x36_5_3_LH2_wf_cycle0_cycle1/log_Iteration_" + std::to_string(i+1) + ".txt";
        // F_value[i] = extractF(filename);
        // cout<<F_value[i]<<endl;
        // h_F_iter->SetBinContent(i+1, F_value[i]);
        // h_F_iter->SetBinError(i+1, 0.00000000001);

        // Matrices
        hh_mat[i] = (TH2D*)infile_debug[i]->Get(Form("hh_mat_%d", i+1));
        hh_mat_inv[i] = (TH2D*)infile_debug[i]->Get(Form("hh_mat_inv_%d", i+1));
        hh_mat_mult[i] = (TH2D*)infile_debug[i]->Get(Form("hh_mat_mult_%d", i+1));
    }

    TCanvas *c_epsilon = new TCanvas("c_epsilon", "", 1000, 1000);
    // h_epsilon_iter->GetYaxis()->SetRangeUser(0, 15);
    format_Mark(h_epsilon_iter, 20, kBlue, 1.5);
    format_Line(h_epsilon_iter, 1, kBlue, 3);
    h_epsilon_iter->Draw("lpe");
    c_epsilon->Print(Form("Result/%s/plots/SummyAndDebugging.pdf)", filename.Data()), "Title:AvgEpsilon_iter");

    // TCanvas *c_F = new TCanvas("c_F", "", 1000, 1000);
    // h_F_iter->GetYaxis()->SetRangeUser(0, 15);
    // format_Mark(h_F_iter, 20, kBlue, 1.5);
    // format_Line(h_F_iter, 1, kBlue, 3);
    // h_F_iter->Draw("lpe");
    // c_F->Print(Form("Result/%s/plots/SummyAndDebugging.pdf)", filename.Data()), "Fvalue_iter");

    TCanvas *c_mat = new TCanvas("c_mat", "", 1000, 1000);
    gStyle->SetPalette(kRainBow);
    for(int i = 0; i<nIter; i++){
        c_mat->cd();
        // c_mat->SetLogz();
        hh_mat[i]->Draw("colz");
        // c_mat->Print(Form("Result/%s/plots/Matrix_iter%d.pdf", filename.Data(), i+1));
        c_mat->Clear();
        hh_mat_inv[i]->Draw("colz");
        // c_mat->Print(Form("Result/%s/plots/Matrix_invert_iter%d.pdf", filename.Data(), i+1));
        c_mat->Clear();
        hh_mat_mult[i]->Draw("colz");
        // c_mat->Print(Form("Result/%s/plots/Matrix_multiply_iter%d.pdf", filename.Data(), i+1));
        c_mat->Clear();
    }

    if(doCyclePlot){ // compare calibration results between different number of cycles
        TGaxis::SetMaxDigits(4);
        TString prefix_2_summary[nTest];
        TString filename_summary[nTest];
        for(int itest = 0; itest < nTest; itest++){
            prefix_2_summary[itest] = Form("cycle%d", startcycle); 
            for(int ic = startcycle+1; ic < startcycle+itest+1; ic++) prefix_2_summary[itest]+=Form("_cycle%d", ic);
            filename_summary[itest] = Form("x%s_%s_%s_%s", Kine.Data(), Tar.Data(), prefix_1.Data(), prefix_2_summary[itest].Data());
        }
        
        TFile *infile_summary_stat[nTest];
        for(int itest = 0; itest < nTest; itest++){
            infile_summary_stat[itest] = TFile::Open(Form("Result/%s/CalibSummary.root", filename_summary[itest].Data()));
        }

        TH1F *h_pi0MeanvsRun_new[nTest];
        TH1F *h_pi0SigmvsRun_new[nTest];
        TH2F *h_coEff_summary[nTest];
        TH2F *hh_coefVariation[nTest];
        for(int itest = 0; itest < nTest; itest++){
            h_pi0MeanvsRun_new[itest] = (TH1F*)infile_summary_stat[itest]->Get("h_pi0MeanvsRun_new");
            h_pi0SigmvsRun_new[itest] = (TH1F*)infile_summary_stat[itest]->Get("h_pi0SigmvsRun_new");
            h_coEff_summary[itest] = (TH2F*)infile_summary_stat[itest]->Get("h_coEff");

            hh_coefVariation[itest] = new TH2F(Form("hh_coefVariation_%d", itest), Form("Variation of coefficients, %d cycles compare to %d cycles (%%);column number;row number", itest+1, itest), 30, -0.5, 29.5, 36, -0.5, 35.5);
        }
        
        // Fill the variation plot
        for(int itest = 1; itest < nTest; itest++){
            for(int iblk = 0; iblk < 1080; iblk++){
                Int_t icol = iblk%30;
                Int_t irow = (iblk-icol)/30;
                Double_t c1 = h_coEff_summary[itest-1]->GetBinContent(iblk+1);
                Double_t c2 = h_coEff_summary[itest]->GetBinContent(iblk+1);
                Double_t vari = abs(100*(c2-c1)/c1);
                hh_coefVariation[itest]->SetBinContent(icol+1, irow+1, vari);
            }
        }

        TLegend *leg = new TLegend(0.55, 0.85-0.05*nTest, 0.85, 0.85);
        for(int itest = 0; itest < nTest; itest++){
            leg->AddEntry(h_pi0MeanvsRun_new[itest], Form("%d cycles combined", itest+1), "lpe");
        }

        Color_t color[7] = {kBlack, kRed, kBlue, kGreen+2, kOrange+3, kMagenta, kViolet+2};
        TCanvas *c = new TCanvas("c", "", 1000, 1000);
        format_Legend(leg);

        // Range of the plots
        Double_t mean_min = 0;
        Double_t mean_max = 0;
        Int_t mean_nbins = h_pi0MeanvsRun_new[0]->GetNbinsX();
        for(int itest = 0; itest < nTest; itest++){
            for(int ibin = 1; ibin<= mean_nbins; ibin++){
                if(mean_min == 0) mean_min = h_pi0MeanvsRun_new[itest]->GetBinContent(ibin);
                if(h_pi0MeanvsRun_new[itest]->GetBinContent(ibin) > 0){
                    if(h_pi0MeanvsRun_new[itest]->GetBinContent(ibin) < mean_min) mean_min = h_pi0MeanvsRun_new[itest]->GetBinContent(ibin);
                    if(h_pi0MeanvsRun_new[itest]->GetBinContent(ibin) > mean_max) mean_max = h_pi0MeanvsRun_new[itest]->GetBinContent(ibin);
                }
            }
        }
        for(int itest = 0; itest < nTest; itest++){ // Draw
            c->cd();
            c->SetGrid();
            
            format_Line(h_pi0MeanvsRun_new[itest], 1, color[itest], 3);
            format_Mark(h_pi0MeanvsRun_new[itest], 20, color[itest], 1.5);
            if(itest == 0){
                h_pi0MeanvsRun_new[itest]->GetYaxis()->SetRangeUser(mean_min-1, mean_max+1);
                h_pi0MeanvsRun_new[itest]->Draw("lpe");
            }
            else h_pi0MeanvsRun_new[itest]->Draw("lpe same");
            leg->Draw();
        }
        c->Print(Form("Result/%s/plots/CompareCycles.pdf(", filename.Data()), "Title:pi0MeanSummary");
        c->Clear();

        // Range of the plots
        Double_t sigm_min = 0;
        Double_t sigm_max = 0;
        Int_t sigm_nbins = h_pi0SigmvsRun_new[0]->GetNbinsX();
        for(int itest = 0; itest < nTest; itest++){
            for(int ibin = 1; ibin<= sigm_nbins; ibin++){
                if(sigm_min == 0) sigm_min = h_pi0SigmvsRun_new[itest]->GetBinContent(ibin);
                if(h_pi0SigmvsRun_new[itest]->GetBinContent(ibin) > 0){
                    if(h_pi0SigmvsRun_new[itest]->GetBinContent(ibin) < sigm_min) sigm_min = h_pi0SigmvsRun_new[itest]->GetBinContent(ibin);
                    if(h_pi0SigmvsRun_new[itest]->GetBinContent(ibin) > 0 
                    && h_pi0SigmvsRun_new[itest]->GetBinContent(ibin) > sigm_max) sigm_max = h_pi0SigmvsRun_new[itest]->GetBinContent(ibin);
                }
            }
        }
        for(int itest = 0; itest < nTest; itest++){ // Draw
            format_Legend(leg);
            c->cd();
            c->SetGrid();

            format_Line(h_pi0SigmvsRun_new[itest], 1, color[itest], 3);
            format_Mark(h_pi0SigmvsRun_new[itest], 20, color[itest], 1.5);
            if(itest == 0){
                h_pi0SigmvsRun_new[itest]->GetYaxis()->SetRangeUser(0.8*sigm_min, 1.2*sigm_max);
                h_pi0SigmvsRun_new[itest]->Draw("lpe");
            }
            else h_pi0SigmvsRun_new[itest]->Draw("lpe same");
            leg->Draw();
        }
        c->Print(Form("Result/%s/plots/CompareCycles.pdf", filename.Data()), "Title:pi0SigmSummary");
        c->Clear();

        for(int itest = 1; itest < nTest; itest++){
            c->cd();
            c->SetCanvasSize(2000, 1000);
            hh_coefVariation[itest]->GetZaxis()->SetRangeUser(0, 5);
            hh_coefVariation[itest]->Draw("colz");
            gStyle->SetPaintTextFormat(".1f");
            hh_coefVariation[itest]->Draw("text same");
            if(itest == nTest-1) c->Print(Form("Result/%s/plots/CompareCycles.pdf)", filename.Data()), Form("Title:coefVariation_%dvs%dcycles", itest+1, itest));
            else c->Print(Form("Result/%s/plots/CompareCycles.pdf", filename.Data()), Form("Title:coefVariation_%dvs%dcycles", itest+1, itest));
            c->Clear();
        } // draw plots
    } // if do summary plots
}