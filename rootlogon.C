{
    cout<<"rootlogon.C has been loaded"<<endl;
    // Load the NPS software path
    //gSystem->Load("libDVCS.so");
    
    char *inc_nps_soft = gSystem->ExpandPathName("/group/nps/hhuang/software/NPS_SOFT");
    gInterpreter->AddIncludePath(inc_nps_soft);
    delete [] inc_nps_soft;
}
