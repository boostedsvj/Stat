void printws(string fname, string space="w"){
	TFile* file = TFile::Open(fname.c_str());
	auto w = (RooWorkspace*)file->Get(space.c_str());
	w->Print();
}
