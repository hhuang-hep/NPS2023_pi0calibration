// #include "/group/nps/hhuang/analysis/MyHeader/Analysis.h"
#include "/group/nps/hhuang/analysis/MyHeader/MyDVCSDB.h"
#include "/group/nps/hhuang/software/NPS_SOFT/TDVCSDB.h"

ofstream logfile;

void MakeList(int kinc_param1, string kinc_param2, int target_flag = -1)
{
    gSystem->Load("/group/nps/hhuang/software/NPS_SOFT/libDVCS.so");
    TDVCSDB *db = new TDVCSDB("dvcs", "clrlpc", 3306, "hhuang", "");

    // Get the run range of each kinematics____________________________________
    vector<int> minrun, maxrun;
    GetRunRange(kinc_param1, kinc_param2, minrun, maxrun);
    TString kinName = Form("x%d_%s", kinc_param1, kinc_param2.c_str());

    Int_t nSubKine = minrun.size();
    if(nSubKine != maxrun.size()){
        cout << "Error: minrun and maxrun vectors have different sizes!" <<endl;
        return;
    }

    // Make a list of production runs_____________________________________________________________
    Int_t Prod_Flag = 0;
    Int_t Quality_Flag = -2;
    Double_t Target_amu = 0;
    Int_t nSegment = -1;
    vector<int> prod_runNb_list, prod_nSeg_list; // list of good production runs should include LH2, LD2, dummy
    vector<int> prod_targType_list; // list of target type of these runs
    for(int i = 0; i < nSubKine; i++){
        for(int run_number = minrun[i]; run_number <= maxrun[i]; run_number++){
            Prod_Flag = *db->GetEntry_i("RUN_flag_IsProdRun", run_number);
            Quality_Flag = *db->GetEntry_i("RUN_flag_Quality", run_number);
            Target_amu = *db->GetEntry_d("TARGET_param_Amu", run_number);
            nSegment = *db->GetEntry_i("DAQ_param_nSegment", run_number);
            
            if(Prod_Flag == 1 && Quality_Flag == 0){
                if(0 < Target_amu && Target_amu < 1.5){ // LH2
                    prod_runNb_list.push_back(run_number);
                    prod_nSeg_list.push_back(nSegment);
                    prod_targType_list.push_back(0);
                }
                if(1.5 < Target_amu && Target_amu < 2.5){ // LD2
                    prod_runNb_list.push_back(run_number);
                    prod_nSeg_list.push_back(nSegment);
                    prod_targType_list.push_back(1);
                }
                if(26.93 < Target_amu && Target_amu < 27.03){ // dummy
                    prod_runNb_list.push_back(run_number);
                    prod_nSeg_list.push_back(nSegment);
                    prod_targType_list.push_back(2);
                }
            }// Select good production runs
        }// loop over run numbers
    }

    // Make a output list for each cycle_______________________________________________

    logfile<<"Start making production list for KinC_"<<kinName.Data()<<endl;

    bool get_LH2 = false;
    bool get_LD2 = false;
    bool get_dummy = false;
    Int_t last_targ_flag = -1;
    vector<vector<int>> v_cycle_runNb; vector<int> v_runNb;
    vector<vector<int>> v_cycle_nSeg; vector<int> v_nSeg;
    vector<vector<int>> v_cycle_targ; vector<int> v_targ;
    for(int i = 0; i < prod_runNb_list.size(); i++){
        if(prod_targType_list.at(i) == 0){
            get_LH2 = true;
            last_targ_flag = prod_targType_list.at(i); // record the last target found in this cycle
            logfile<<"Found LH2 run "<<prod_runNb_list.at(i)<<" with "<<prod_nSeg_list.at(i)<<" segments."<<endl;

            v_runNb.push_back(prod_runNb_list.at(i)); // record the run number for this cycle
            v_nSeg.push_back(prod_nSeg_list.at(i)); // record the number of segment for this cycle
            v_targ.push_back(prod_targType_list.at(i)); // record the target type for this cycle
        }
        if(prod_targType_list.at(i) == 1){
            get_LD2 = true;
            last_targ_flag = prod_targType_list.at(i); // record the last target found in this cycle
            logfile<<"Found LD2 run "<<prod_runNb_list.at(i)<<" with "<<prod_nSeg_list.at(i)<<" segments."<<endl;

            v_runNb.push_back(prod_runNb_list.at(i)); // record the run number for this cycle
            v_nSeg.push_back(prod_nSeg_list.at(i)); // record the number of segment for this cycle
            v_targ.push_back(prod_targType_list.at(i)); // record the target type for this cycle
        }
        if(prod_targType_list.at(i) == 2){
            get_dummy = true;
            last_targ_flag = prod_targType_list.at(i); // record the last target found in this cycle
            logfile<<"Found dummy run "<<prod_runNb_list.at(i)<<" with "<<prod_nSeg_list.at(i)<<" segments."<<endl;

            v_runNb.push_back(prod_runNb_list.at(i)); // record the run number for this cycle
            v_nSeg.push_back(prod_nSeg_list.at(i)); // record the number of segment for this cycle
            v_targ.push_back(prod_targType_list.at(i)); // record the target type for this cycle
        }

        // if found all 3 target and about to change cycle or finish the whole list
        bool isLastRun = (i + 1 == prod_runNb_list.size());
        bool targetWillChange = !isLastRun && (prod_targType_list.at(i+1) != last_targ_flag);
        
        if(get_LH2 && get_LD2 && get_dummy && (isLastRun || targetWillChange)){

            // put the current into the vector of cycles
            v_cycle_runNb.push_back(v_runNb); // run list
            v_cycle_nSeg.push_back(v_nSeg); // nseg list
            v_cycle_targ.push_back(v_targ); // target list

            // clean up for the next cycle
            v_runNb.clear();
            v_nSeg.clear();
            v_targ.clear();
            get_LH2 = false;
            get_LD2 = false;
            get_dummy = false;
        }
    }

    // After looping all runs, put leftover runs into a new cycle
    if(!v_runNb.empty()){
        logfile << "Leftover runs found (incomplete cycle), adding as extra cycle." << endl;
        logfile << "This cycle contains targets: ";

        bool hasLH2 = false, hasLD2 = false, hasDummy = false;
        for(auto t : v_targ){
            if(t == 0) hasLH2 = true;
            if(t == 1) hasLD2 = true;
            if(t == 2) hasDummy = true;
        }

        logfile << "["
                << (hasLH2   ? "LH2 "   : "")
                << (hasLD2   ? "LD2 "   : "")
                << (hasDummy? "Dummy " : "")
                << "]" << endl;

        logfile << "Runs in this incomplete cycle:" << endl;
        for(int j = 0; j < v_runNb.size(); j++){
            logfile << "  Run " << v_runNb[j] << " (Target: "
                    << (v_targ[j] == 0 ? "LH2" : (v_targ[j] == 1 ? "LD2" : "Dummy"))
                    << ", Segments: " << v_nSeg[j] << ")" << endl;
        }
        
        // push this incomplete cycle
        v_cycle_runNb.push_back(v_runNb);
        v_cycle_nSeg.push_back(v_nSeg);
        v_cycle_targ.push_back(v_targ);

        // cleanup (optional, since function ends)
        v_runNb.clear();
        v_nSeg.clear();
        v_targ.clear();
    }

    logfile<<"============================================================"<<endl;

    // Loop over the vector of cycles to output the list
    get_LH2 = false;
    get_LD2 = false;
    get_dummy = false;
    ofstream fProdList_LH2, fProdList_LD2, fProdList_dummy;
    for(int i = 0; i < v_cycle_runNb.size(); i++){
        for(int j = 0; j < v_cycle_runNb[i].size(); j++){
            Int_t run = v_cycle_runNb[i][j];
            Int_t nseg = v_cycle_nSeg[i][j];
            Int_t targ = v_cycle_targ[i][j];
            if(targ == 0){
                if(!get_LH2){
                    // Make the list file when get the first run of that target
                    system(Form("rm -f ListPerCycle/%s_%d_cycle%d.txt", kinName.Data(), targ, i));
                    fProdList_LH2.open(Form("ListPerCycle/%s_%d_cycle%d.txt", kinName.Data(), targ, i));
                    if(!(fProdList_LH2.is_open())) cout<<"Error: can't open the list for LH2"<<endl;
                    get_LH2 = true;
                }
                fProdList_LH2<<run<<" "<<nseg<<endl;
            } // targ == LH2
            if(targ == 1){
                if(!get_LD2){
                    // Make the list file when get the first run of that target
                    system(Form("rm -f ListPerCycle/%s_%d_cycle%d.txt", kinName.Data(), targ, i));
                    fProdList_LD2.open(Form("ListPerCycle/%s_%d_cycle%d.txt", kinName.Data(), targ, i));
                    if(!(fProdList_LD2.is_open())) cout<<"Error: can't open the list for LD2"<<endl;
                    get_LD2 = true;
                }
                fProdList_LD2<<run<<" "<<nseg<<endl;
            } // targ == LD2
            if(targ == 2){
                if(!get_dummy){
                    // Make the list file when get the first run of that target
                    system(Form("rm -f ListPerCycle/%s_%d_cycle%d.txt", kinName.Data(), targ, i));
                    fProdList_dummy.open(Form("ListPerCycle/%s_%d_cycle%d.txt", kinName.Data(), targ, i));
                    if(!(fProdList_dummy.is_open())) cout<<"Error: can't open the list for dummy"<<endl;
                    get_dummy = true;
                }
                fProdList_dummy<<run<<" "<<nseg<<endl;
            } // targ == dummy
        } // loop over the runs in the cycle

        // close output files of this cycle
        fProdList_LH2.close();
        fProdList_LD2.close();
        fProdList_dummy.close();
        // reset for next cycle
        get_LH2 = false;
        get_LD2 = false;
        get_dummy = false;
    } // cycle loop
}

void MakeListPerCycle()
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

    system("rm -f ListPerCycle/log_cycleList.txt");
    logfile.open("ListPerCycle/log_cycleList.txt");
    
    for(int i=0; i<54; i++){
        MakeList(array_par1[i], array_par2[i]);
    }
    cout<<"Finished making all calibration lists"<<endl;
}