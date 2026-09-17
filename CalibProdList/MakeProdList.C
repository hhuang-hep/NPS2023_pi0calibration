// #include "/group/nps/hhuang/analysis/MyHeader/Analysis.h"
#include "/group/nps/hhuang/analysis/MyHeader/MyDVCSDB.h"
#include "/group/nps/hhuang/software/NPS_SOFT/TDVCSDB.h"

void MakeList(int kinc_param1, string kinc_param2, int target_flag = -1)
{
    gSystem->Load("/group/nps/hhuang/software/NPS_SOFT/libDVCS.so");
    TDVCSDB *db = new TDVCSDB("dvcs", "clrlpc", 3306, "hhuang", "");

    vector<int> minrun, maxrun;
    GetRunRange(kinc_param1, kinc_param2, minrun, maxrun);
    TString kinName = Form("x%d_%s", kinc_param1, kinc_param2.c_str());

    Int_t nSubKine = minrun.size();
    if(nSubKine != maxrun.size()){
        cout << "Error: minrun and maxrun vectors have different sizes!" << endl;
        return;
    }

    // Specify the target type?
    bool isLH2, isLD2, isDummy;
    if(target_flag == -1){
        isLH2 = true;
        isLD2 = true;
        isDummy = true;
    }
    else if(target_flag == 0){
        isLH2 = true;
        isLD2 = false;
        isDummy = false;
    }
    else if(target_flag == 1){
        isLH2 = false;
        isLD2 = true;
        isDummy  = false;
    }
    else if(target_flag == 2){
        isLH2 = false;
        isLD2 = false;
        isDummy  = true;
    }
    else{
        cout << "Error: Invalid target flag!" << endl;
     
        return;
    }
    
    if(isLH2 && !isLD2 && !isDummy){
        // target_flag = 0;
        cout<<"Start making production list for KinC_"<<kinName.Data()<<endl; 
        for(int i = 0; i < nSubKine; i++){
            cout<<"Min Run: "<<minrun[i]<<", Max Run: "<<maxrun[i]<<endl;
        }
        cout<<"Target type: LH2"<< endl;
    }
    else if(!isLH2 && isLD2 && !isDummy){
        cout<<"Start making production list for KinC_"<<kinName.Data()<<endl; 
        for(int i = 0; i < nSubKine; i++){
            cout<<"Min Run: "<<minrun[i]<<", Max Run: "<<maxrun[i]<<endl;
        }
        cout<<"Target type: LD2"<< endl;
    }
    else if(!isLH2 && !isLD2 && isDummy){
        cout<<"Start making production list for KinC_"<<kinName.Data()<<endl; 
        for(int i = 0; i < nSubKine; i++){
            cout<<"Min Run: "<<minrun[i]<<", Max Run: "<<maxrun[i]<<endl;
        }
        cout<<"Target type: Dummy"<< endl;
    }
    else if(isLH2 && isLD2 && isDummy){
        cout<<"Start making production list for KinC_"<<kinName.Data()<<endl; 
        for(int i = 0; i < nSubKine; i++){
            cout<<"Min Run: "<<minrun[i]<<", Max Run: "<<maxrun[i]<<endl;
        }
        cout<<"Target type: LH2+LD2+Dummy"<< endl;
    }
    else{
        cout << "Error: Invalid target type!" << endl;
        return;
    }

    system(Form("rm -f %s_%d_ProdList.txt", kinName.Data(), target_flag));
    ofstream fProdList(Form("%s_%d_ProdList.txt", kinName.Data(), target_flag));

    // This list is for simulation and data analysis,
    // which only inludes "good" "production" runs, so DO NOT change flags below
    Int_t Prod_Flag = 0;
    Int_t Quality_Flag = -2;
    Double_t Target_amu = 0;
    Int_t nSegment = -1;
    for(int i = 0; i < nSubKine; i++){
        for(int run_number = minrun[i]; run_number <= maxrun[i]; run_number++){
            Prod_Flag = *db->GetEntry_i("RUN_flag_IsProdRun", run_number);
            Quality_Flag = *db->GetEntry_i("RUN_flag_Quality", run_number);
            Target_amu = *db->GetEntry_d("TARGET_param_Amu", run_number);
            nSegment = *db->GetEntry_i("DAQ_param_nSegment", run_number);
            
            if(Prod_Flag == 1 && Quality_Flag == 0){ // to select good production runs
                if(target_flag == -1){ // LH2+LD2
                    if(0 < Target_amu && Target_amu < 2.5) fProdList << run_number <<" "<< nSegment << endl; 
                }
                else if(target_flag == 0){// LH2
                    if(0 < Target_amu && Target_amu < 1.5) fProdList << run_number <<" "<< nSegment << endl; 
                }
                else if(target_flag == 1){// LD2
                    if(1.5 < Target_amu && Target_amu < 2.5) fProdList << run_number <<" "<< nSegment << endl;
                }
                else if(target_flag == 2){// Dummy
                    if(26.93 < Target_amu && Target_amu < 27.03) fProdList << run_number <<" "<< nSegment << endl;
                }
            }
        }
    }
}

void MakeProdList()
{
    // 54 Kinematics
    const Int_t nKine = 54;
    Int_t array_par1[nKine] = {
        // KinC_x25 (5)
        25, 25, 25, 25, 25,
        // KinC_x36 (16)
        36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36,
        // KinC_x50 (19)
        50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50,
        // KinC_x60 (14)
        60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60
    };
    string array_par2[nKine] = {
        // KinC_x25
        "1_1", "1_2", "3", "4_1", "4_2",
        // KinC_x36
        "1", "2_1", "2_2", "2p", "2pp", "3_1", "3_2", "4", "5_1", "5_2", "5_3", "5_4", "5p", "6_1", "6_2", "6_3",
        // KinC_x50
        "0_1", "0_2", "0a_1", "0a_2", "0b_1", "0b_2", "1_1", "1_2", "1p", "2_1", "2_2", "2_3", "2p", "2pp", "3_1", "3_2", "3_3", "3p", "3pp",
        // KinC_x60
        "1", "2b", "2_1", "2_2", "2p", "3_1", "3_2", "3_3", "3p", "3a", "3b", "4a_1", "4a_2", "4b"
    };

    for(int i=0; i<54; i++){
        MakeList(array_par1[i], array_par2[i], -1);
        MakeList(array_par1[i], array_par2[i], 0);
        MakeList(array_par1[i], array_par2[i], 1);
        MakeList(array_par1[i], array_par2[i], 2);
    }
    cout<<"Finished making all production lists"<<endl;
}