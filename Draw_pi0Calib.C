#include "/group/nps/hhuang/analysis/MyHeader/MyDB.h"
#include "/group/nps/hhuang/analysis/MyHeader/DrawSetting.h"
#include "/group/nps/hhuang/analysis/MyHeader/Analysis.h"

void Draw_pi0Calib(TString Kine, int TargetFlag)
{
    Global_setting();

    const Int_t nIter = 8;

    TString Tar;
    if(TargetFlag == 0) Tar = "LH2";
    else if(TargetFlag == 1) Tar = "LD2";
    else if(TargetFlag == -1) Tar = "LH2_LD2";
    else{
        cout<<"ERROR: Unknown target!!!"<<endl;
        return;
    }

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

    // Connect to the database and get kinematics variables________________________________________________
    gSystem->Load("/group/nps/hhuang/software/NPS_SOFT/libDVCS.so");
    TDVCSDB *db = new TDVCSDB("dvcs", "clrlpc", 3306, "hhuang", "");

    // Elastic coefficients
    Double_t *coefElas = new Double_t[1080]; // Get the elastic coefficients in NPS numbering Scheme
    coefElas = db->GetEntry_d("CALO_calib_ElasCoef", first_run);

    // Get pi0 correction factors
    TString filename = Form("%s_%s", Kine.Data(), Tar.Data());
    ifstream fcorr_old(Form("Result/%s_pass2_v2/corr_pi0.txt", filename.Data())); // In simulation numbering scheme

    Double_t corr_pi0[1080];
    for(int iblk = 0; iblk < 1080; iblk++) fcorr_old>>corr_pi0[bnConv_OldToNew(iblk)]; // convert to NPS numbering scheme

    // Histograms
    TH1F *h_coEff = new TH1F("h_coEff", "Calibration coefficients after correction;PMT number;Corrected coefficients", 1080, -0.5,  1079.5);
    h_coEff->SetTitle("Calibration coefficients;PMT number;coefficients");

    TH2F *hh_coEff = new TH2F("hh_coEff", "Calibration coefficients; Column number; Row number; ", 30, -0.5, 29.5, 36, -0.5, 35.5);

    TH2F *hh_iblk = new TH2F("hh_iblk", "", 30, -0.5, 29.5, 36, -0.5, 35.5);

    // system(Form("rm -f coef_final/for_hcana/coef_pi0Calib_%s.txt", filename.Data()));
    // ofstream fout_coef_hcana(Form("coef_final/for_hcana/coef_pi0Calib_%s.txt", filename.Data()));
    // fout_coef_hcana<<"nps_cal_arr_gain_cor = ";

    system(Form("rm -f Result/%s_pass2_v2/coef_pi0Calib_temp.txt", filename.Data()));
    ofstream fout_coef(Form("Result/%s_pass2_v2/coef_pi0Calib_temp.txt", filename.Data()));
    for(int iblk = 0; iblk < 1080; iblk++){
        Int_t icol = iblk%30;
        Int_t irow = (iblk-icol)/30;
        hh_iblk->SetBinContent(icol+1, irow+1, iblk);

        Double_t coef_final = coefElas[iblk]*corr_pi0[iblk];

        // Fill the calibration coefficients into histogram
        h_coEff->SetBinContent(iblk+1, coef_final);
        hh_coEff->SetBinContent(icol+1, irow+1, coef_final);

        // Output coefficients
        fout_coef<<coef_final<<endl;
        // if((iblk+1) % 30 == 0) fout_coef_hcana<<coef_final<<", \n           ";
        // else fout_coef_hcana<<coef_final<<", ";
    }

    // Coefficients for each run----->Need reconstruction again for further correction
    // for(int irun = 0; irun < nRun; irun++){
    //     system((Form("cp coef_pi0Calib_temp.txt coef_final/coef_pi0Calib_%d.txt", runList[irun])));
    // }

    // system("rm -f coef_pi0Calib_temp.txt");

    TFile *infile[nIter+1];
    TH1 *h_pi0M;
    TF1 *f_fit;
    TLegend *leg_pi0M;
    TCanvas *c_pi0M;
    TH1F *h_mass_iter = new TH1F("h_mass_iter", "M_{#pi^{0}} vs. iteration;Iteration;M_{#pi^{0}} [MeV/c^{2}]",nIter+1, -0.5, nIter+0.5);
    TH1F *h_width_iter = new TH1F("h_width_iter", "#sigma_{#pi^{0}} vs. iteration;Iteration;#sigma_{#pi^{0}} [MeV/c^{2}]",nIter+1, -0.5, nIter+0.5);
    
    for(int iIter = 0; iIter < nIter+1; iIter++){
        infile[iIter] = TFile::Open(Form("Result/%s_pass2_v2/f_Mgg_%d.root", filename.Data(), iIter));
        h_pi0M = (TH1F*)infile[iIter]->Get(Form("h_pi0M_%d", iIter));
        h_pi0M->Rebin(2);
        Double_t meanfit = h_pi0M->GetBinCenter(h_pi0M->GetMaximumBin());
        Double_t minfit = meanfit - 0.008;
        Double_t maxfit = meanfit + 0.008;
        f_fit = new TF1(Form("f_fit_%d", iIter), "gausn(0)", minfit, maxfit);
        // f_fit = new TF1(Form("f_fit_%d", iIter), "gausn(0)+pol1(3)", 0.08, 0.2);

        f_fit->SetParLimits(0, 1, 10000);
        f_fit->SetParLimits(1, 0.115, 0.15);
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

        leg_pi0M = new TLegend(0.40, 0.65, 0.85, 0.90);
        leg_pi0M->AddEntry(h_pi0M, "Data", "lpe");
        leg_pi0M->AddEntry(f_fit, "Fit (gausn+pol1)", "l");
        leg_pi0M->AddEntry((TObject*)0, Form("#chi^{2}/NDF = %.2f", chi2/ndf), "");
        leg_pi0M->AddEntry((TObject*)0, Form("N_{#pi^{0}} = %.0f #pm %.0f", N, N_err), "");
        leg_pi0M->AddEntry((TObject*)0, Form("M_{#pi^{0}} = %.2f #pm %.2f MeV", mean, mean_err), "");
        leg_pi0M->AddEntry((TObject*)0, Form("#sigma_{#pi^{0}} = %.2f #pm %.2f MeV", sigm, sigm_err), "");

        c_pi0M = new TCanvas(Form("c_pi0M_%d", iIter), "", 1000, 1000);
        h_pi0M->GetYaxis()->SetRangeUser(0, 1.8*h_pi0M->GetMaximum());
        format_Mark(h_pi0M, 20, kBlack, 1.5);
        format_Line(h_pi0M, 1, kBlack, 3);
        h_pi0M->Draw("lpe");

        // Draw the pi0 mass before calibration
        if(iIter == nIter){
            TH1F *h_pi0M_0 = (TH1F*)infile[0]->Get("h_pi0M_0");
            // h_pi0M_0->Rebin(2);

            format_Line(h_pi0M_0, 1, kGreen+2, 3);
            h_pi0M_0->Draw("same hist");
        }

        format_Line(f_fit, 1, kRed, 5);
        f_fit->Draw("C same");

        format_Legend(leg_pi0M);
        leg_pi0M->Draw();

        c_pi0M->SaveAs(Form("Result/%s_pass2_v2/pi0M_iter_%d.png", filename.Data(), iIter));
    }

    TLine *l_mass = new TLine(-0.5, 1000*0.1349766, nIter+0.5, 1000*0.1349766);

    TLegend *leg_mass = new TLegend(0.45, 0.75, 0.85, 0.85);
    leg_mass->SetHeader("#bf{E_{photons} > 1.4 GeV}");
    leg_mass->AddEntry(l_mass, "M_{#pi^{0}} = 0.1349766 GeV/c^{2}", "l");

    TLegend *leg_width = new TLegend(0.45, 0.80, 0.85, 0.85);
    leg_width->SetHeader("#bf{E_{photons} > 1.4 GeV}");

    TCanvas *c_mass = new TCanvas("c_mass", "", 1000, 1000);
    h_mass_iter->GetYaxis()->SetRangeUser(100, 150);
    format_Mark(h_mass_iter, 20, kBlue, 1.5);
    format_Line(h_mass_iter, 1, kBlue, 3);
    h_mass_iter->Draw("lpe");

    format_Line(l_mass, 7, kRed, 3);
    l_mass->Draw("l same");

    format_Legend(leg_mass);
    leg_mass->Draw();

    c_mass->SaveAs(Form("Result/%s_pass2_v2/Mass_iter.png", filename.Data()));

    TCanvas *c_width = new TCanvas("c_width", "", 1000, 1000);
    h_width_iter->GetYaxis()->SetRangeUser(0.8*h_width_iter->GetMinimum(), 1.2*h_width_iter->GetMaximum());
    format_Mark(h_width_iter, 20, kBlue, 1.5);
    format_Line(h_width_iter, 1, kBlue, 3);
    h_width_iter->Draw("lpe");

    format_Legend(leg_width);
    leg_width->Draw();

    c_width->SaveAs(Form("Result/%s_pass2_v2/Width_iter.png", filename.Data()));

    TCanvas *c_coEff = new TCanvas("c_coEff", "", 1000, 1000);
    format_Line(h_coEff, 1, kBlue, 3);
    h_coEff->Draw("hist");
    c_coEff->SaveAs(Form("Result/%s_pass2_v2/coefficient.png", filename.Data()));
    
    h_coEff->GetYaxis()->SetRangeUser(0, 2);
    c_coEff->SaveAs(Form("Result/%s_pass2_v2/coefficient_zoomin.png", filename.Data()));

    TCanvas *c_coEff2D = new TCanvas("c_coEff2D", "", 2000, 1000);
    hh_coEff->Draw("colz");
    hh_iblk->Draw("text same");
    c_coEff2D->SaveAs(Form("Result/%s_pass2_v2/coefficient2D.png", filename.Data()));

    hh_coEff->GetZaxis()->SetRangeUser(0, 2);
    hh_iblk->Draw("text same");
    c_coEff2D->SaveAs(Form("Result/%s_pass2_v2/coefficient2D_zoomin.png", filename.Data()));
}