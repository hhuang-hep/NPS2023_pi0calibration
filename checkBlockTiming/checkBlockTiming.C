#include "/group/nps/hhuang/analysis/MyHeader/DrawSetting.h"
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TList.h>
#include <TString.h>
#include <iostream>

void checkBlockTiming(TString kine = "x36_4"){
    Global_setting();
    
    TString calibTreeDir = Form("/group/nps/hhuang/analysis/DVCS_NPS2023/DVCS_analysis/pi0Calib_wf/calibTree/%s", kine.Data());
    TSystemDirectory dir("calibTree", calibTreeDir.Data());
    TList* files = dir.GetListOfFiles();

    if (!files) {
        cerr << "Cannot open directory " << calibTreeDir.Data() <<endl;
        return;
    }

    TIter next(files);
    TSystemFile* file;

    TFile *infile;
    TH1F *h_temp;
    TH1F *h_NpsTime[1080];
    for(int iblk = 0; iblk < 1080; iblk++){
        h_NpsTime[iblk] = new TH1F(Form("h_NpsTime_blk%d", iblk), "NPS timing with offset correction (ampl. > 10 mV);Pulse time [ns];Number of events", 20000, -100, 100);
    }

    while ((file = (TSystemFile*)next())) {
        TString fname = file->GetName();
        cout << fname <<endl;

        // skip . and ..
        if (fname == "." || fname == "..") continue;

        // skip sub-directories
        if (file->IsDirectory()) continue;

        infile = TFile::Open(Form("%s/%s", calibTreeDir.Data(), fname.Data()));
        for(int i = 0; i < 1080; i++){
            h_temp = (TH1F*)infile->Get(Form("h_NpsTime_blk%d", i));
            h_NpsTime[i]->Add(h_temp);
        }

        infile->Close();
    }

    TF1 *f_fitNpsTime[1080];
    for(int iblk = 0; iblk < 1080; iblk++) f_fitNpsTime[iblk] = new TF1(Form("f_fitNpsTime_blk%d", iblk), "gausn(0)+[3]*sin(TMath::Pi()*x[0]+[4])+[5]", -10, 10);

    for(int iblk = 0; iblk < 1080; iblk++){
        f_fitNpsTime[iblk]->SetParLimits(0, 1, 10000);
        f_fitNpsTime[iblk]->SetParameter(1, h_NpsTime[iblk]->GetBinCenter(h_NpsTime[iblk]->GetMaximumBin()));
        f_fitNpsTime[iblk]->SetParLimits(2, 0.00001, 1);

        h_NpsTime[iblk]->Rebin(20);
        h_NpsTime[iblk]->Fit(f_fitNpsTime[iblk], "RN");

        
    }
    
    TH1F *h_fitResult = new TH1F("h_fitResult", "NPS timing fit results;Block number;Timing peak (ns)", 1080, -0.5, 1079.5);
    for(int iblk = 0; iblk < 1080; iblk++){
        h_fitResult->SetBinContent(iblk+1, f_fitNpsTime[iblk]->GetParameter(1));
        h_fitResult->SetBinError(iblk+1, f_fitNpsTime[iblk]->GetParError(1));
    }

    TCanvas *c = new TCanvas("c", "", 1000, 1000);
    c->SetGridx();
    for(int iblk = 0; iblk < 1080; iblk++){
        h_NpsTime[iblk]->GetXaxis()->SetRangeUser(-10, 10);
        format_Line(h_NpsTime[iblk], 1, kBlack, 3);
        format_Mark(h_NpsTime[iblk], 20, kBlack, 1.5);
        h_NpsTime[iblk]->Draw("lpe");
        format_Line(f_fitNpsTime[iblk], 1, kRed, 3);
        f_fitNpsTime[iblk]->Draw("same C");
        if(iblk == 0) c->Print(Form("plots/%s_wfTime.pdf(", kine.Data()), Form("Title:h_NpsTime_blk%d", iblk));
        else c->Print(Form("plots/%s_wfTime.pdf", kine.Data()), Form("Title:h_NpsTime_blk%d", iblk));
        c->Clear();
    }

    c->SetCanvasSize(2000,1000);
    c->SetGridx(0);
    c->SetGridy();
    format_Line(h_fitResult, 1, kBlue, 3);
    format_Mark(h_fitResult, 20, kBlue, 1.5);
    h_fitResult->Draw("lpe");
    c->Print(Form("plots/%s_wfTime.pdf", kine.Data()), "Title:NpsTimeFitResult");
    h_fitResult->GetYaxis()->SetRangeUser(-10, 10);
    c->Print(Form("plots/%s_wfTime.pdf", kine.Data()), "Title:NpsTimeFitResult_zoomin_1");
    h_fitResult->GetYaxis()->SetRangeUser(-0.5, 0.5);
    c->Print(Form("plots/%s_wfTime.pdf)", kine.Data()), "Title:NpsTimeFitResult_zoomin_2");
}