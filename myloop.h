#include <vector>

using namespace std;

#define MUON_MASS    0.10565837
#define PION_MASS    0.13957018
#define KAON_MASS    0.493677
#define KSHORT_MASS  0.497614
#define KSTAR_MASS   0.89594
#define PHI_MASS     1.019455
#define JPSI_MASS    3.096916
#define PSI2S_MASS   3.686109
#define PROTON_MASS  0.938272046
#define LAMBDA_MASS  1.115683
#define BP_MASS      5.27926
#define B0_MASS      5.27958
#define BS_MASS      5.36677
#define BC_MASS      6.2756
#define LAMBDAB_MASS 5.6195

// HLT paths to be booked
enum {
    HLT_Dimuon0_Jpsi_Muon_v1,            // 0
    HLT_Dimuon10_Jpsi_Barrel_v1,         // 1
    HLT_Dimuon16_Jpsi_v1,                // 2
    HLT_Dimuon20_Jpsi_v1,                // 3
    HLT_DoubleMu4_3_Jpsi_Displaced_v1,   // 4
    HLT_DoubleMu4_JpsiTrk_Displaced_v2,  // 5
    N_HLT_BOOKINGS
};

const char HLT_paths[N_HLT_BOOKINGS][64] = {
    "HLT_Dimuon0_Jpsi_Muon_v1",          // 0
    "HLT_Dimuon10_Jpsi_Barrel_v1",       // 1
    "HLT_Dimuon16_Jpsi_v1",              // 2
    "HLT_Dimuon20_Jpsi_v1",              // 3
    "HLT_DoubleMu4_3_Jpsi_Displaced_v1", // 4
    "HLT_DoubleMu4_JpsiTrk_Displaced_v2" // 5
};

class ReducedBranches{
    public:
    UInt_t     runNum;
    UInt_t     evtNum;
    //
    UInt_t nMCB0;
    vector<float> *MCmupPx, *MCmupPy, *MCmupPz ;
    vector<float> *MCmumPx, *MCmumPy, *MCmumPz ;
    vector<float> *MCpionPx, *MCpionPy, *MCpionPz ;
    vector<float> *MCkaonPx, *MCkaonPy, *MCkaonPz ;
    /*
    vector<float> *MCmupPx = 0; vector<float> *MCmupPy = 0; vector<float> *MCmupPz = 0;
    vector<float> *MCmumPx = 0; vector<float> *MCmumPy = 0; vector<float> *MCmumPz = 0;
    vector<float> *MCpionPx = 0; vector<float> *MCpionPy = 0; vector<float> *MCpionPz = 0;
    vector<float> *MCkaonPx = 0; vector<float> *MCkaonPy = 0; vector<float> *MCkaonPz = 0;
    */
    
    TBranch        *b_MCmupPx, *b_MCmupPy, *b_MCmupPz ;
    TBranch        *b_MCmumPx, *b_MCmumPy, *b_MCmumPz ;
    TBranch        *b_MCpionPx, *b_MCpionPy, *b_MCpionPz ;
    TBranch        *b_MCkaonPx, *b_MCkaonPy, *b_MCkaonPz ;
    /*
    TBranch *b_MCmupPx = 0; TBranch *b_MCmupPy = 0; TBranch *b_MCmupPz = 0;
    TBranch *b_MCmumPx = 0; TBranch *b_MCmumPy = 0; TBranch *b_MCmumPz = 0;
    TBranch *b_MCpionPx = 0; TBranch *b_MCpionPy = 0; TBranch *b_MCpionPz = 0;
    TBranch *b_MCkaonPx = 0; TBranch *b_MCkaonPy = 0; TBranch *b_MCkaonPz = 0;
    */
    void regTree(TTree *root){
      root->Branch("runNum",&runNum,"runNum/i");
      root->Branch("evtNum",&evtNum,"evtNum/i");
      //
      /*
      MCmupPx = 0; MCmupPy = 0; MCmupPz = 0;
      MCmumPx = 0; MCmumPy = 0; MCmumPz = 0;
      MCpionPx = 0; MCpionPy = 0; MCpionPz = 0;
      MCkaonPx = 0; MCkaonPy = 0; MCkaonPz = 0;
      */
      // 
      root->Branch("nMCB0", &nMCB0);
      root->Branch("MCmupPx", &MCmupPx);
      root->Branch("MCmupPy", &MCmupPy);
      root->Branch("MCmupPz", &MCmupPz);
      root->Branch("MCmumPx", &MCmumPx);
      root->Branch("MCmumPy", &MCmumPy);
      root->Branch("MCmumPz", &MCmumPz);
      root->Branch("MCpionPx", &MCpionPx);
      root->Branch("MCpionPy", &MCpionPy);
      root->Branch("MCpionPz", &MCpionPz);
      root->Branch("MCkaonPx", &MCkaonPx);
      root->Branch("MCkaonPy", &MCkaonPy);
      root->Branch("MCkaonPz", &MCkaonPz);
    }

    void setbranchadd(TTree *root){      
      root->SetBranchAddress("runNum",&runNum);
      root->SetBranchAddress("evtNum",&evtNum);
      //
      cout <<"\nbefore setting vectors branches" <<endl;
      root->SetBranchAddress("nMCB0", &nMCB0);
      root->SetBranchAddress("MCmupPx", &MCmupPx);
      root->SetBranchAddress("MCmupPy", &MCmupPy);
      root->SetBranchAddress("MCmupPz", &MCmupPz);
      root->SetBranchAddress("MCmumPx", &MCmumPx);
      root->SetBranchAddress("MCmumPy", &MCmumPy);
      root->SetBranchAddress("MCmumPz", &MCmumPz);
      root->SetBranchAddress("MCpionPx", &MCpionPx);
      root->SetBranchAddress("MCpionPy", &MCpionPy);
      root->SetBranchAddress("MCpionPz", &MCpionPz);
      root->SetBranchAddress("MCkaonPx", &MCkaonPx);
      root->SetBranchAddress("MCkaonPy", &MCkaonPy);
      root->SetBranchAddress("MCkaonPz", &MCkaonPz);
      /*
      root->SetBranchAddress("MCmupPx", &MCmupPx, &b_MCmupPx);
      root->SetBranchAddress("MCmupPy", &MCmupPy, &b_MCmupPy);
      root->SetBranchAddress("MCmupPz", &MCmupPz, &b_MCmupPz);
      root->SetBranchAddress("MCmumPx", &MCmumPx, &b_MCmumPx);
      root->SetBranchAddress("MCmumPy", &MCmumPy, &b_MCmumPy);
      root->SetBranchAddress("MCmumPz", &MCmumPz, &b_MCmumPz);
      root->SetBranchAddress("MCpionPx", &MCpionPx, &b_MCpionPx);
      root->SetBranchAddress("MCpionPy", &MCpionPy, &b_MCpionPy);
      root->SetBranchAddress("MCpionPz", &MCpionPz, &b_MCpionPz);
      root->SetBranchAddress("MCkaonPx", &MCkaonPx, &b_MCkaonPx);
      root->SetBranchAddress("MCkaonPy", &MCkaonPy, &b_MCkaonPy);
      root->SetBranchAddress("MCkaonPz", &MCkaonPz, &b_MCkaonPz);
      */
      cout <<"\nafter setting vectors branches" <<endl;
    }
};
