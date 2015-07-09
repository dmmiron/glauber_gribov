void MakeJets() {
    gROOT->LoadMacro("density.cxx");
    gROOT->LoadMacro("sampling.cxx");
    gROOT->ProcessLine("SweepJets(2000, 0, 0, 16, 1, 0, 91, 30, \"jetpairs_2k/\", true, -10, -10, 10, 10)");
    //gROOT->ProcessLine("SweepJets(20000, 0, 8, 16, 1, 0, 91, 100, \"jets/\", false, -10, -10, 10, 10)");
}

