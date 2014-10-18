




#define NUMPHIS 3
#define NUMETAS 13
#define NUMDEPTHS 2

struct TCalibLedInfo
{
	int numChs;
	int iphi[50];
	int ieta[50];
	int cBoxChannel[50];
	vector<string> *cBoxString;
	int nTS[50];
	double pulse[50][10];
};

struct THFInfo
{
	int numChs;
	int numTS;
	int iphi[2000];
	int ieta[2000];
	int depth[2000];
	double pulse[2000][50];
};

bool isGarbage(double pulse[50]);

struct H2Triggers
{
	int ped;
	int led;
	int laser;
	int beam;
};

struct BCData
{
	double s1adc;
	double s2adc;
	double s3adc;
	double s4adc;
};

struct H2MapChannel
{
	int i;
	int crate;
	int slot;
	char topBottom;
	int dcc;
	int spigot;
	int fiber;
	int fiberCh;
	string subDet;
	string subSubDet;
	int ieta;
	int iphi;
	int depth;
	int channelId;
};

//
//	Was lazy, init map HF(iphi, ieta, depth) -> channelId
//
int HF2ChId[NUMPHIS][NUMETAS][NUMDEPTHS];
HF2ChId[1][0][0] = 1;
HF2ChId[1][1][0] = 2;
HF2ChId[1][2][0] = 3;
HF2ChId[1][3][0] = 4;
HF2ChId[1][4][0] = 5;
HF2ChId[1][5][0] = 6;
HF2ChId[1][6][0] = 7;
HF2ChId[1][7][0] = 8;
HF2ChId[1][8][0] = 9;
HF2ChId[1][9][0] = 10;
HF2ChId[1][10][0] = 11;
HF2ChId[1][0][1] = 12;
HF2ChId[1][1][1] = 13;
HF2ChId[1][2][1] = 14;
HF2ChId[1][3][1] = 15;
HF2ChId[1][4][1] = 16;
HF2ChId[1][5][1] = 17;
HF2ChId[1][6][1] = 18;
HF2ChId[1][7][1] = 19;
HF2ChId[1][8][1] = 20;
HF2ChId[1][9][1] = 21;
HF2ChId[1][10][1] = 22;
HF2ChId[2][0][0] = 23;
HF2ChId[2][1][0] = 24;



vector<double> *wcX[5];
vector<double> *wcY[5];

void genH2Plots(int verb, string inFileName, string outFileName)
{
	TCalibLedInfo calibInfo;
	THFInfo hfInfo;
	H2Triggers triggers;
	BCData bc;

	ofstream log_beam("log.beam");
	ofstream log_ped("log.ped");
	ofstream log_beam_strange("log.beam.strange");
	ofstream log_beam_excl("log.beam.excl");

	TFile *in = new TFile(inFileName.c_str());

	in->cd("HFData");
	TTree *treeHF = (TTree*)gDirectory->Get("Events");
	treeHF->SetBranchAddress("numChs", &hfInfo.numChs);
	treeHF->SetBranchAddress("numTS", &hfInfo.numTS);
	treeHF->SetBranchAddress("iphi", &hfInfo.iphi);
	treeHF->SetBranchAddress("ieta", &hfInfo.ieta);
	treeHF->SetBranchAddress("depth", &hfInfo.depth);
	treeHF->SetBranchAddress("pulse", hfInfo.pulse);

	in->cd("Triggers");
	TTree *treeTriggers = (TTree*)gDirectory->Get("Events");
	treeTriggers->SetBranchAddress("ped", &triggers.ped);
	treeTriggers->SetBranchAddress("led", &triggers.led);
	treeTriggers->SetBranchAddress("beam", &triggers.beam);

	in->cd("WCData");
	TTree *treeWC = (TTree*)gDirectory->Get("Events");
	treeWC->SetBranchAddress("xA", &(wcX[0]));
	treeWC->SetBranchAddress("yA", &(wcY[0]));
	treeWC->SetBranchAddress("xB", &(wcX[1]));
	treeWC->SetBranchAddress("yB", &(wcY[1]));
	treeWC->SetBranchAddress("xC", &(wcX[2]));
	treeWC->SetBranchAddress("yC", &(wcY[2]));
	treeWC->SetBranchAddress("xD", &(wcX[3]));
	treeWC->SetBranchAddress("yD", &(wcY[3]));
	treeWC->SetBranchAddress("xE", &(wcX[4]));
	treeWC->SetBranchAddress("yE", &(wcY[4]));

	in->cd("BeamCounters");
	TTree *treeBC = (TTree*)gDirectory->Get("Events");
	treeBC->SetBranchAddress("s1adc", &bc.s1adc);
	treeBC->SetBranchAddress("s2adc", &bc.s2adc);
	treeBC->SetBranchAddress("s3adc", &bc.s3adc);
	treeBC->SetBranchAddress("s4adc", &bc.s4adc);


	TFile *out = new TFile(outFileName.c_str(), "recreate");
	out->mkdir("PedEvents");
	out->mkdir("DataEvents");

	out->cd("PedEvents");
	gDirectory->mkdir("Histos");
	gDirectory->mkdir("AvgPulses");

	out->cd("DataEvents");
	gDirectory->mkdir("Histos");
	gDirectory->mkdir("AvgPulses");
	gDirectory->mkdir("XY");

	TH1D *hSignal[2][NUMPHIS][NUMETAS][NUMDEPTHS];
	TProfile *pSignalPulses[2][NUMPHIS][NUMETAS][NUMDEPTHS];
	TH2D *hXY[5];
	TH2D *hXY_wCut[5];
	char wcNames[5] = {'A', 'B', 'C', 'D', 'E'};
	char nameXY[200];

	out->cd("DataEvents/XY");
	for (int i=0; i<5; i++)
	{
		sprintf(nameXY, "XY%c", wcNames[i]);
		hXY[i] = new TH2D(nameXY, nameXY, 200, -100, 100, 200, -100, 100);
		sprintf(nameXY, "XY%c_wCut", wcNames[i]);
		hXY_wCut[i] = new TH2D(nameXY, nameXY, 200, -100, 100, 200, -100, 100);
	}

	char histName[200];
	for (int itype=0; itype<2; itype++)
	{
		if (itype==0)
			out->cd("PedEvents");
		else if (itype==1)
			out->cd("DataEvents");
		for (int iiphi=0; iiphi<NUMPHIS; iiphi++)
		{
			for (int iieta=0; iieta<NUMETAS; iieta++)
			{
				for (int idepth=0; idepth<NUMDEPTHS; idepth++)
				{
					if (HF2ChId[iiphi][iieta][idepth]==0)
					{
						gDirectory->cd("Histos");
						sprintf(histName, "Signal_IPHI%d_IETA%d_D%d",
								2*iiphi+1, iieta+29, idepth+1);
						hSignal[itype][iiphi][iieta][idepth] = new TH1D(
								histName, histName, 5000, 0, 1000);
		
						gDirectory->cd("../AvgPulses");
						sprintf(histName, "SignalAvgPulse_IPHI%d_IETA%d_D%d",
								2*iiphi+1, iieta+29, idepth+1);
						pSignalPulses[itype][iiphi][iieta][idepth] = new TProfile(
								histName, histName, 30, 0, 30);

						gDirectory->cd("..");
	
						continue;
					}


					gDirectory->cd("Histos");
					sprintf(histName, "Signal_Ch%d",
							HF2ChId[iiphi][iieta][idepth]);
					hSignal[itype][iiphi][iieta][idepth] = new TH1D(
							histName, histName, 5000, 0, 1000);
	
					gDirectory->cd("../AvgPulses");
					sprintf(histName, "SignalAvgPulse_Ch%d",
							HF2ChId[iiphi][iieta][idepth]);
					pSignalPulses[itype][iiphi][iieta][idepth] = new TProfile(
							histName, histName, 30, 0, 30);

					gDirectory->cd("..");
				}
			}
		}
	}

	cout << "### Starting the Analysis..." << endl;
	int numEvents = treeHF->GetEntries();
	cout << "### Total Number of Events=" << numEvents << endl;
	for (int iEvent=0; iEvent<numEvents; iEvent++)
	{
		treeHF->GetEntry(iEvent);
		treeTriggers->GetEntry(iEvent);
		treeWC->GetEntry(iEvent);
		treeBC->GetEntry(iEvent);

//		cout << "HEREHERE" << endl;

		//
		//	Fill Wire Chambers XY Profiles
		//
/*		for (int i=0; i<5; i++)
		{
			double xx = wcX[i]->at(0);
			double yy = wcY[i]->at(0);
			hXY[i]->Fill(xx, yy);
		}
*/

		if (iEvent%1000==0)
			cout << "### Event=" << iEvent << endl;

		if (triggers.ped==1 && (triggers.beam==0 || triggers.led==0))
		{
			
			//
			//	Ped Trigger Only -> Get the Pedestal 
			//	itype=0
			//
			for (int iCh=0; iCh<hfInfo.numChs; iCh++)
			{
				int iphi = hfInfo.iphi[iCh]; int iiphi=(iphi-1)/2;
				int ieta = hfInfo.ieta[iCh]; int iieta=ieta-29;
				int depth = hfInfo.depth[iCh]; int idepth = depth-1;
				int iiphi = (iphi-1)/2;
				int iieta = ieta-29;
				int idepth = depth-1;
	
				Double_t totSigSum = 0;
				if (verb>1 && iEvent<=10)
					cout << "### Digi->Pulse: " << iphi << "  " << ieta << "  "
						<< depth << endl;

				log_ped << HF2ChId[iiphi][iieta][idepth] << endl;
				for (int iTS=0; iTS<hfInfo.numTS; iTS++)
				{
					totSigSum += hfInfo.pulse[iCh][iTS];
					pSignalPulses[0][iiphi][iieta][idepth]->Fill(iTS, 
							hfInfo.pulse[iCh][iTS]);
	
					log_ped << hfInfo.pulse[iCh][iTS] << "  ";
				}
				log_ped << endl;

				hSignal[0][iiphi][iieta][idepth]->Fill(totSigSum);
			}
			continue;
		}


		//
		//	For Beam or LED Events
		//	itype=1 for that
		//
		for (int i=0; i<5; i++)
		{
			if (wcX[i]->size()==0 || wcY[i]->size()==0)
				continue;
			double xx = wcX[i]->at(0);
			double yy = wcY[i]->at(0);

			hXY[i]->Fill(xx, yy);

			if (bc.s2adc>60)
				hXY_wCut[i]->Fill(xx, yy);
		}

		for (int iCh=0; iCh<hfInfo.numChs; iCh++)
		{
			int iphi = hfInfo.iphi[iCh]; int iiphi = (iphi-1)/2;
			int ieta = hfInfo.ieta[iCh]; int iieta = ieta-29;
			int depth = hfInfo.depth[iCh]; int idepth = depth-1;
			int iiphi = (iphi-1)/2;
			int iieta = ieta-29;
			int idepth = depth-1;

			if (isGarbage(hfInfo.pulse[iCh]))
			{
				log_beam_excl << HF2ChId[iiphi][iieta][idepth] << endl;
				for (int i=0; i<10; i++)
					log_beam_excl << hfInfo.pulse[iCh][i] << "  ";
				log_beam_excl << endl;
				continue;
			}

			Double_t totSigSum = 0;
			if (verb>1 && iEvent<=10)
				cout << "### Digi->Pulse: " << iphi << "  " << ieta << "  "
					<< depth << endl;
			log_beam << HF2ChId[iiphi][iieta][idepth] << endl;
			for (int iTS=0; iTS<hfInfo.numTS; iTS++)
			{
				totSigSum += hfInfo.pulse[iCh][iTS];
				pSignalPulses[1][iiphi][iieta][idepth]->Fill(iTS, 
						hfInfo.pulse[iCh][iTS]);

				log_beam << hfInfo.pulse[iCh][iTS] << "  ";
			}
			log_beam << endl;
			hSignal[1][iiphi][iieta][idepth]->Fill(totSigSum);

			if (totSigSum>100 && totSigSum<200)
			{
				log_beam_strange << HF2ChId[iiphi][iieta][idepth] << endl;
				for (int i=0; i<10; i++)
					log_beam_strange << hfInfo.pulse[iCh][i] << "  ";
				log_beam_strange << endl;
			}
			
		}
	}

	out->Write();
	out->Close();

	return;
}

bool isGarbage(double pulse[50])
{
	double sum = 0;
	for (int i=0; i<3; i++)
		sum += pulse[i];

	sum /= 3;
	if (sum>50)
		return true;
	else 
		return false;

	return false;
}
















