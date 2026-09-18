Int_t GetKinematics(int run, double &HMS_mom, double &Target_amu, double &HMS_ang, double &SHMS_ang, double &Beam_energy)
{
    ifstream inputFile;
    inputFile.open(Form("/work/hallc/nps/nps-ana/REPORT_OUTPUT_pass1/COIN/SKIM/skim_NPS_HMS_report_%d_-1.report", run));

    if(inputFile.is_open()){
        cout<<"Get the report file of run"<<run<<endl;
    }
    if(!inputFile.is_open()){
        cout<<"Error: Report file of run"<<run<<" is not found!!"<<endl;
        return 0;
    }
    

    string line;
    while(std::getline(inputFile, line))
    {   
        if(line.find("Momentum") != std::string::npos) {
            std::sscanf(line.c_str(), "Momentum : %lf", &HMS_mom);
        }
        else if(line.find("Target AMU") != std::string::npos){
            std::sscanf(line.c_str(), "Target AMU  : %lf", &Target_amu);
        }
        else if(line.find("SHMS Theta") != std::string::npos){
            std::sscanf(line.c_str(), "SHMS Theta  : %lf", &SHMS_ang);
        }
        else if(line.find("HMS Theta") != std::string::npos){
            std::sscanf(line.c_str(), "HMS Theta  : %lf", &HMS_ang);
        }
        else if(line.find("Beam Energy") != std::string::npos) {
            std::sscanf(line.c_str(), "Beam Energy : %lf", &Beam_energy);
            break;
        }
    }
    inputFile.close();
    return 1;
}

Int_t GetPsFactor(int run, int &trig, int &ps)
{
    ifstream inputFile;
    inputFile.open(Form("/work/hallc/nps/nps-ana/REPORT_OUTPUT_pass1/COIN/SKIM/skim_NPS_HMS_report_%d_-1.report", run));

    if(inputFile.is_open()){
        cout<<"Get the report file of run"<<run<<endl;
    }
    if(!inputFile.is_open()){
        cout<<"Error: Report file of run"<<run<<" is not found!!"<<endl;
        return 0;
    }

    string line;
    while(std::getline(inputFile, line))
    {   
        if(line.find("Ps1_factor") != std::string::npos) {
            trig = 1;
            std::sscanf(line.c_str(), "Ps1_factor = %d", &ps);
            if(ps != -1) break;
        }

        if(line.find("Ps2_factor") != std::string::npos) {
            trig = 2;
            std::sscanf(line.c_str(), "Ps2_factor = %d", &ps);
            if(ps != -1) break;
        }

        if(line.find("Ps3_factor") != std::string::npos) {
            trig = 3;
            std::sscanf(line.c_str(), "Ps3_factor = %d", &ps);
            if(ps != -1) break;
        }

        if(line.find("Ps4_factor") != std::string::npos) {
            trig = 4;
            std::sscanf(line.c_str(), "Ps4_factor = %d", &ps);
            if(ps != -1) break;
        }

        if(line.find("Ps5_factor") != std::string::npos) {
            trig = 5;
            std::sscanf(line.c_str(), "Ps5_factor = %d", &ps);
            if(ps != -1) break;
        }

        if(line.find("Ps6_factor") != std::string::npos) {
            trig = 6;
            std::sscanf(line.c_str(), "Ps6_factor = %d", &ps);
            if(ps != -1) break;
        }
    }

    inputFile.close();
    return 1;
}

Int_t GetCharge(int run, double &charge)
{
    ifstream inputFile;
    inputFile.open(Form("/work/hallc/nps/nps-ana/REPORT_OUTPUT_pass1/COIN/SKIM/skim_NPS_HMS_report_%d_-1.report", run));

    if(inputFile.is_open()){
        cout<<"Get the report file of run"<<run<<endl;
    }
    if(!inputFile.is_open()){
        cout<<"Error: Report file of run"<<run<<" is not found!!"<<endl;
        return 0;
    }

    string line;
    while(std::getline(inputFile, line))
    {   
        if(line.find("BCM4A Beam Cut Charge") != std::string::npos) {
            std::sscanf(line.c_str(), "BCM4A Beam Cut Charge: %lf", &charge);
        }
    }

    inputFile.close();
    return 1;
}

Int_t GetTrigger(int run, int &Ps1_factor, int &Ps2_factor, int &Ps3_factor, int &Ps4_factor, int &Ps5_factor, int &Ps6_factor)
{
    ifstream inputFile;
    inputFile.open(Form("/work/hallc/nps/nps-ana/REPORT_OUTPUT_pass1/COIN/SKIM/skim_NPS_HMS_report_%d_-1.report", run));
    
    if(inputFile.is_open()){
        cout<<"Get the report file of run"<<run<<endl;
    }
    if(!inputFile.is_open()){
        cout<<"Error: Report file of run"<<run<<" is not found!!"<<endl;
        return 0;
    }

    string line;
    while(std::getline(inputFile, line))
    {   
        if (line.find("Ps1_factor") != std::string::npos) {
            cout<<line<<endl;
            std::sscanf(line.c_str(), "Ps1_factor = %d", &Ps1_factor);
        } 
        else if (line.find("Ps2_factor") != std::string::npos) {
            cout<<line<<endl;
            std::sscanf(line.c_str(), "Ps2_factor = %d", &Ps2_factor);
        } 
        else if (line.find("Ps3_factor") != std::string::npos) {
            cout<<line<<endl;
            std::sscanf(line.c_str(), "Ps3_factor = %d", &Ps3_factor);
        } 
        else if (line.find("Ps4_factor") != std::string::npos) {
            cout<<line<<endl;
            std::sscanf(line.c_str(), "Ps4_factor = %d", &Ps4_factor);
        } 
        else if (line.find("Ps5_factor") != std::string::npos) {
            cout<<line<<endl;
            std::sscanf(line.c_str(), "Ps5_factor = %d", &Ps5_factor);
        } 
        else if (line.find("Ps6_factor") != std::string::npos) {
            cout<<line<<endl;
            std::sscanf(line.c_str(), "Ps6_factor = %d", &Ps6_factor);
            break;
        }
    }
    inputFile.close();
    return 1;
}

TString GetKineName(double HMS_mom, double HMS_angle, double SHMS_angle, double Beam_energy, double NPS_dist)
{
    TString KinC_main, KinC_sub;

    if(HMS_mom == 6.667){
        KinC_main = "50-2";
        if(SHMS_angle == 36.88) KinC_sub = "";
        else if(SHMS_angle == 35.446) KinC_sub = "p";
        else if(SHMS_angle == 38.31) KinC_sub = "pp";
    }
    else if(HMS_mom == 4.042){
        KinC_main = "36_2";
        if(SHMS_angle == 30.664) KinC_sub = "";
        else if(SHMS_angle == 28.76) KinC_sub = "p";
        else if(SHMS_angle == 32.9) KinC_sub = "pp";
    }
    else if(HMS_mom == 6.117){
        KinC_main = "36_3";
        KinC_sub = "";
    }
    else if(HMS_mom == 4.637){
        KinC_main = "36_5";
        if(SHMS_angle == 28.417) KinC_sub = "";
        else if(SHMS_angle == 30.30) KinC_sub = "p";
    }
    else if(HMS_mom == 5.878){
        KinC_main = "60_3";
        if(SHMS_angle == 35.02) KinC_sub = "";
        else if(SHMS_angle == 34.018) KinC_sub = "p";
        else if(SHMS_angle == 33.015) KinC_sub = "a";
        else if(SHMS_angle == 36.45) KinC_sub = "b";
    }
    else if(HMS_mom == 5.038){
        KinC_main = "60_4";
        if(SHMS_angle == 30.375) KinC_sub = "a";
        else if(SHMS_angle == 33.812) KinC_sub = "b";
    }
    else if(HMS_mom == 2.416){
        KinC_main = "36_6";
        KinC_sub = "";
    }
    else if(HMS_mom == 4.149){
        KinC_main = "25_4";
        KinC_sub = "";
    }
    else if(HMS_mom == 5.253 && NPS_dist == 300){
        KinC_main = "50_4";
        KinC_sub = "";
    }
    else if(HMS_mom == 5.253 && NPS_dist == 400){
        KinC_main = "50_3";
        if(SHMS_angle == 29.742) KinC_sub = "p";
        else if(SHMS_angle == 33.179) KinC_sub = "pp";
    }
    else if(HMS_mom == 4.726){
        KinC_main = "50_1";
        if(SHMS_angle == 35.284 || SHMS_angle == 35.285) KinC_sub = "";
        else if(SHMS_angle == 33.38 || SHMS_angle == 33.379) KinC_sub = "p";
    }
    else if(HMS_mom == 3.803){
        KinC_main = "60_2";
        if(SHMS_angle == 32.867 || SHMS_angle == 32.869) KinC_sub = "";
        else if(SHMS_angle == 28.75 || SHMS_angle == 28.76) KinC_sub = "p";
    }

    TString KinC_name = Form("KinC_x%s%s", KinC_main.Data(), KinC_sub.Data());
    return KinC_name;
}

TString GetTarName(double Target_amu)
{   
    TString name;
    if(abs(Target_amu-1) < 0.5) name = "LH2";
    if(abs(Target_amu-2) < 0.5) name = "LD2";
    if(abs(Target_amu-26.98) < 0.5) name = "Dummy";

    return name;
}

TString GetFile_ElasCoef(int runNb)
{
    TString coefDir = "/work/hallc/nps/hhuang/elasCalib/coef";
    TString fileName;
    if(1561 <= runNb && runNb <= 1968) fileName = "elasCoef_npsSchm_3";
    else if(1983 <= runNb && runNb <= 2854) fileName = "elasCoef_npsSchm_4";
    else if(2921 <= runNb && runNb <= 3727) fileName = "elasCoef_npsSchm_7";
    else if(3728 <= runNb && runNb <= 3882) fileName = "elasCoef_npsSchm_8";
    else if(3899 <= runNb) fileName = "elasCoef_npsSchm_9";
    else{
        TString error = "Error: No file was found!!!!";
        cout<<"Please check if there are coefficients for this run"<<endl;
        return error;
    }

    TString file = Form("%s/%s.txt", coefDir.Data(), fileName.Data());
    return file;
}

TString GetListName_pi0Calib(int runNb)
{
    TString fileName = "";
    
    if(3728 <= runNb && runNb <= 3737)  fileName = "36-5_1_0";

    else if((3990 <= runNb && runNb <= 3998) || (4017 <= runNb && runNb <= 4023)) fileName = "36-5p_1_0";
    else if((4002 <= runNb && runNb <= 4003) || runNb == 4009)  fileName = "36-5p_2_0";
    else if(4010 <= runNb && runNb <= 4016) fileName = "36-5p_2_1";
    else if(4036 <= runNb && runNb <= 4038) fileName = "36-5p_2_2";
    else if(4039 <= runNb && runNb <= 4045) fileName = "36-5p_2_3";

    else if(1662 <= runNb && runNb <= 1675) fileName = "50-2_2_0";

    else if(4224 <= runNb && runNb <= 4233) fileName = "60-4a_1_0";
    else if(4271 <= runNb && runNb <= 4278) fileName = "60-4a_1_1";
    else if(4325 <= runNb && runNb <= 4334) fileName = "60-4a_1_2";
    else if(4373 <= runNb && runNb <= 4380) fileName = "60-4a_1_3";
    
    else if((4242 <= runNb && runNb <= 4252)) fileName = "60-4a_2_0";
    else if((4283 <= runNb && runNb <= 4299)) fileName = "60-4a_2_1";
    else if((4337 <= runNb && runNb <= 4348)) fileName = "60-4a_2_2";
    else if((4423 <= runNb && runNb <= 4482)) fileName = "60-4a_2_3";
    else if((4510 <= runNb && runNb <= 4550)) fileName = "60-4a_2_4";

    else{
        cout<<"Error: No file was found!!!! Please check if there is a list or coefficient file for this run"<<endl;
    }

    return fileName;
}

TString GetFile_Pi0Coef(int runNb)
{   
    TString coefDir = ("/work/hallc/nps/hhuang/pi0Calib/coef_final");

    TString fileName = GetListName_pi0Calib(runNb);
    TString file = Form("%s/%s.txt", coefDir.Data(), fileName.Data());

    if(fileName != "") return file;
    else return fileName;
}

// Double_t GetBeamEnergy(int runNb)
// {
//     Double_t eBeam; // beam energy [GeV]
//     if(1572 <= runNb && runNb <= 1752) eBeam = 10.538;
//     if(3728 <= runNb && runNb < 9999) eBeam = 10.538;

//     return eBeam;
// }

// Double_t GetCaloAngle(int runNb)
// {
//     Double_t ShmsAngle; // degree
//     Double_t caloOffset = 16.3; // degree
//     if(1572 <= runNb && runNb <= 1752) ShmsAngle = 36.880; // 50-2
//     if(3728 <= runNb && runNb <= 3782) ShmsAngle = 28.497; // 36-5
//     if(3783 <= runNb && runNb <= 3806) ShmsAngle = 35.441; // 50-2'
//     if(3807 <= runNb && runNb <= 3832) ShmsAngle = 37.995; // 50-2"
//     if(3833 <= runNb && runNb <= 3857) ShmsAngle = 35.438; // 50-2'
//     if(3858 <= runNb && runNb <= 3882) ShmsAngle = 37.995; // 50-2"
//     if(3990 <= runNb && runNb <= 4045) ShmsAngle = 30.305; // 36-5'
//     if(4047 <= runNb && runNb <= 4075) ShmsAngle = 33.015; // 60-3a
//     if(4076 <= runNb && runNb <= 4091) ShmsAngle = 36.445; // 60-3b
//     if(4224 <= runNb && runNb <= 4252) ShmsAngle = 30.373; // 60-4a
//     if(4271 <= runNb && runNb <= 4299) ShmsAngle = 30.370; // 60-4a
//     if(4325 <= runNb && runNb <= 4348) ShmsAngle = 30.378; // 60-4a
//     if(4373 <= runNb && runNb <= 4394) ShmsAngle = 30.372; // 60-4a
//     if(4423 <= runNb && runNb <= 4482) ShmsAngle = 30.373; // 60-4a
//     if(4510 <= runNb && runNb <= 4517) ShmsAngle = 30.375; // 60-4a
//     Double_t caloAngle = (ShmsAngle-caloOffset)*TMath::DegToRad();
//     return caloAngle;
// }