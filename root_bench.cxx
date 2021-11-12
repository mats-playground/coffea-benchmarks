// Compile with:
// $CXX -O3 $(root-config --cflags --libs) -o root_bench root_bench.cxx
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "Math/Vector4D.h"

double query1(const char * filename) {
    ROOT::RDataFrame df("Events", filename);
    auto h = df.Histo1D<float>({"", ";MET (GeV);N_{Events}", 100, 0, 200}, "MET_pt");

    return h->Integral();
}

double query2(const char * filename) {
    ROOT::RDataFrame df("Events", filename);
    auto h = df.Histo1D<ROOT::RVec<float>>({"", ";Jet p_{T} (GeV);N_{Events}", 100, 15, 60}, "Jet_pt");

    return h->Integral();
}

double query3(const char * filename) {
    ROOT::RDataFrame df("Events", filename);
    auto goodJetPt = [](const ROOT::RVec<float> &pt, const ROOT::RVec<float> &eta) { return pt[abs(eta) < 1.0]; };
    auto h = df.Define("goodJet_pt", goodJetPt, {"Jet_pt", "Jet_eta"})
               .Histo1D<ROOT::RVec<float>>({"", ";Jet p_{T} (GeV);N_{Events}", 100, 15, 60}, "goodJet_pt");

    return h->Integral();
}

double query4(const char * filename) {
    ROOT::RDataFrame df("Events", filename);
    auto filter = [](const ROOT::RVec<float> & pt, const ROOT::RVec<float> & eta) {
            return Sum(pt > 40) > 1;
    };
    auto h = df.Filter(filter, {"Jet_pt", "Jet_eta"}, "More than one jet with pt > 40")
               .Histo1D<float>({"", ";MET (GeV);N_{Events}", 100, 0, 200}, "MET_pt");

    return h->Integral();
}

template <typename T> using Vec = const ROOT::RVec<T>&;
using FourVector = ROOT::Math::PtEtaPhiMVector;

auto compute_dimuon_masses(Vec<float> pt, Vec<float> eta, Vec<float> phi, Vec<float> mass, Vec<int> charge)
{
    ROOT::RVec<float> masses;
    // std::cout << pt.size() << std::endl;
    const auto c = ROOT::VecOps::Combinations(pt, 2);
    for (auto i = 0; i < c[0].size(); i++) {
        const auto i1 = c[0][i];
        const auto i2 = c[1][i];
        if (charge[i1] == charge[i2]) continue;
        const FourVector p1(pt[i1], eta[i1], phi[i1], mass[i1]);
        const FourVector p2(pt[i2], eta[i2], phi[i2], mass[i2]);
        masses.push_back((p1 + p2).mass());
    }
    return masses;
};

double query5(const char * filename) {
    ROOT::RDataFrame df("Events", filename);
    auto h = df.Filter([](unsigned int n) { return n >= 2; }, {"nMuon"}, "At least two muons")
               .Define("Dimuon_mass", compute_dimuon_masses,
                       {"Muon_pt", "Muon_eta", "Muon_phi", "Muon_mass", "Muon_charge"})
               .Filter([](const ROOT::RVec<float> &mass) { return Sum(mass > 60 && mass < 120) > 0; }, {"Dimuon_mass"},
                       "At least one dimuon system with mass in range [60, 120]")
               .Histo1D<float>({"", ";MET (GeV);N_{Events}", 100, 0, 200}, "MET_pt");

    return h->Integral();
}

ROOT::RVec<std::size_t> find_trijet(Vec<float> pt, Vec<float> eta, Vec<float> phi, Vec<float> mass)
{
    const auto c = ROOT::VecOps::Combinations(pt, 3);
    const auto make_p4 = [&](std::size_t idx) {
        return ROOT::Math::PtEtaPhiMVector(pt[idx], eta[idx], phi[idx], mass[idx]);
    };

    float trijet_mass = -1;
    float distance = 1e9;
    const auto top_mass = 172.5;
    std::size_t idx = 0;
    for (auto i = 0; i < c[0].size(); i++) {
        auto p1 = make_p4(c[0][i]);
        auto p2 = make_p4(c[1][i]);
        auto p3 = make_p4(c[2][i]);
        const auto tmp_mass = (p1 + p2 + p3).mass();
        const auto tmp_distance = std::abs(tmp_mass - top_mass);
        if (tmp_distance < distance) {
            distance = tmp_distance;
            trijet_mass = tmp_mass;
            idx = i;
        }

    }

    return {c[0][idx], c[1][idx], c[2][idx]};
}

float trijet_pt(Vec<float> pt, Vec<float> eta, Vec<float> phi, Vec<float> mass, Vec<std::size_t> idx)
{
    auto p1 = ROOT::Math::PtEtaPhiMVector(pt[idx[0]], eta[idx[0]], phi[idx[0]], mass[idx[0]]);
    auto p2 = ROOT::Math::PtEtaPhiMVector(pt[idx[1]], eta[idx[1]], phi[idx[1]], mass[idx[1]]);
    auto p3 = ROOT::Math::PtEtaPhiMVector(pt[idx[2]], eta[idx[2]], phi[idx[2]], mass[idx[2]]);
    return (p1 + p2 + p3).pt();
}

double query6(const char * filename) {
    ROOT::RDataFrame df("Events", filename);
    auto df2 = df.Filter([](unsigned int n) { return n >= 3; }, {"nJet"}, "At least three jets")
                 .Define("Trijet_idx", find_trijet, {"Jet_pt", "Jet_eta", "Jet_phi", "Jet_mass"});
    auto h1 = df2.Define("Trijet_pt", trijet_pt, {"Jet_pt", "Jet_eta", "Jet_phi", "Jet_mass", "Trijet_idx"})
                 .Histo1D<float>({"", ";Trijet pt (GeV);N_{Events}", 100, 15, 40}, "Trijet_pt");
    auto h2 = df2.Define("Trijet_leadingBtag",
                         [](const ROOT::RVec<float> &btag, const ROOT::RVec<std::size_t> &idx) { return Max(Take(btag, idx)); },
                         {"Jet_btag", "Trijet_idx"})
                 .Histo1D<float>({"", ";Trijet leading b-tag;N_{Events}", 100, 0, 1}, "Trijet_leadingBtag");

    return h1->Integral() + h2->Integral();
}

ROOT::RVec<int> find_isolated_jets(Vec<float> eta1, Vec<float> phi1, Vec<float> pt2, Vec<float> eta2, Vec<float> phi2)
{
    ROOT::RVec<int> mask(eta1.size(), 1);
    if (eta2.size() == 0) {
        return mask;
    }

    const auto ptcut = pt2 > 10;
    const auto eta2_ptcut = eta2[ptcut];
    const auto phi2_ptcut = phi2[ptcut];
    if (eta2_ptcut.size() == 0) {
        return mask;
    }

    const auto c = ROOT::VecOps::Combinations(eta1, eta2_ptcut);
    for (auto i = 0; i < c[0].size(); i++) {
        const auto i1 = c[0][i];
        const auto i2 = c[1][i];
        const auto dr = ROOT::VecOps::DeltaR(eta1[i1], eta2_ptcut[i2], phi1[i1], phi2_ptcut[i2]);
        if (dr < 0.4) mask[i1] = 0;
    }
    return mask;
}

double query7(const char * filename) {
    ROOT::RDataFrame df("Events", filename);
    auto h = df.Filter([](unsigned int n) { return n > 0; }, {"nJet"}, "At least one jet")
               .Define("goodJet_ptcut", [](const ROOT::RVec<float>& pt) { return pt > 30; }, {"Jet_pt"})
               .Define("goodJet_antiMuon",
                       find_isolated_jets,
                       {"Jet_eta", "Jet_phi", "Muon_pt", "Muon_eta", "Muon_phi"})
               .Define("goodJet_antiElectron",
                       find_isolated_jets,
                       {"Jet_eta", "Jet_phi", "Electron_pt", "Electron_eta", "Electron_phi"})
               .Define("goodJet",
                       [](const ROOT::RVec<int> &pt, const ROOT::RVec<int> &muon, const ROOT::RVec<int> &electron) {
                            return pt && muon && electron;
                       },
                       {"goodJet_ptcut", "goodJet_antiMuon", "goodJet_antiElectron"})
               .Filter([](const ROOT::RVec<int> &good) { return Sum(good) > 0; }, {"goodJet"})
               .Define("goodJet_sumPt",
                       [](const ROOT::RVec<int> &good, const ROOT::RVec<float> &pt) { return Sum(pt[good]); },
                       {"goodJet", "Jet_pt"})
               .Histo1D<float>({"", ";Jet p_{T} sum (GeV);N_{Events}", 100, 15, 200}, "goodJet_sumPt");

    return h->Integral();
}

unsigned int additional_lepton_idx(Vec<float> pt, Vec<float> eta, Vec<float> phi, Vec<float> mass, Vec<int> charge, Vec<int> flavour)
{
    const auto c = Combinations(pt, 2);
    unsigned int lep_idx = -999;
    float best_mass = 99999;
    int best_i1 = -1;
    int best_i2 = -1;
    const auto z_mass = 91.2;
    const auto make_p4 = [&](std::size_t idx) {
        return ROOT::Math::PtEtaPhiMVector(pt[idx], eta[idx], phi[idx], mass[idx]);
    };

    for (auto i = 0; i < c[0].size(); i++) {
        const auto i1 = c[0][i];
        const auto i2 = c[1][i];
        if (charge[i1] == charge[i2]) continue;
        if (flavour[i1] != flavour[i2]) continue;
        const auto p1 = make_p4(i1);
        const auto p2 = make_p4(i2);
        const auto tmp_mass = (p1 + p2).mass();
        if (std::abs(tmp_mass - z_mass) < std::abs(best_mass - z_mass)) {
            best_mass = tmp_mass;
            best_i1 = i1;
            best_i2 = i2;
        }
    }

    if (best_i1 == -1) return lep_idx;

    float max_pt = -999;
    for (auto i = 0; i < pt.size(); i++) {
        if (i != best_i1 && i != best_i2 && pt[i] > max_pt) {
            max_pt = pt[i];
            lep_idx = i;
        }
    }

    return lep_idx;
}

double query8(const char * filename) {
    ROOT::RDataFrame df("Events", filename);
    auto concatF = [](const ROOT::RVec<float> &a, const ROOT::RVec<float> &b) { return Concatenate(a, b); };
    auto concatI = [](const ROOT::RVec<int> &a, const ROOT::RVec<int> &b) { return Concatenate(a, b); };
    auto transverseMass = [](const ROOT::RVec<float> &Lepton_pt, const ROOT::RVec<float> &Lepton_phi,
                             float MET_pt, float MET_phi, unsigned int idx) {
        return sqrt(2.0 * Lepton_pt[idx] * MET_pt * (1.0 - cos(ROOT::VecOps::DeltaPhi(MET_phi, Lepton_phi[idx]))));
    };
    auto h = df.Filter([](unsigned int nElectron, unsigned int nMuon) { return nElectron + nMuon > 2; },
                       {"nElectron", "nMuon"}, "At least three leptons")
               .Define("Lepton_pt", concatF, {"Muon_pt", "Electron_pt"})
               .Define("Lepton_eta", concatF, {"Muon_eta", "Electron_eta"})
               .Define("Lepton_phi", concatF, {"Muon_phi", "Electron_phi"})
               .Define("Lepton_mass", concatF, {"Muon_mass", "Electron_mass"})
               .Define("Lepton_charge", concatI, {"Muon_charge", "Electron_charge"})
               .Define("Lepton_flavour",[](unsigned int nMuon, unsigned int nElectron) {
                                            return Concatenate(ROOT::RVec<int>(nMuon, 0), ROOT::RVec<int>(nElectron, 1));
                                        },
                       {"nMuon", "nElectron"})
               .Define("AdditionalLepton_idx", additional_lepton_idx,
                       {"Lepton_pt", "Lepton_eta", "Lepton_phi", "Lepton_mass", "Lepton_charge", "Lepton_flavour"})
               .Filter([](unsigned int idx) { return idx != -999; }, {"AdditionalLepton_idx"}, "No valid lepton pair found.")
               .Define("TransverseMass", transverseMass,
                       {"Lepton_pt", "Lepton_phi", "MET_pt", "MET_phi", "AdditionalLepton_idx"})

               .Histo1D<double>({"", ";Transverse mass (GeV);N_{Events}", 100, 0, 200}, "TransverseMass");

    return h->Integral();
}

int main(int argc, const char ** argv) {
    if ( argc != 3 ) {
      std::cout << "Usage: PROG ncores filename" << std::endl;
      return 1;
    }
    int ncores = std::stoi(argv[1]);
    const char * filename = argv[2];

    if ( ncores > 1 ) {
      ROOT::EnableImplicitMT(ncores);
    }

    std::cout << "query,ncores,walltime" << std::endl;

    auto tic = std::chrono::steady_clock::now();
    query1(filename);
    auto toc = std::chrono::steady_clock::now();
    std::chrono::duration<double> diff = toc - tic;
    std::cout << "Q1," << ncores << "," << diff.count() << std::endl;

    tic = std::chrono::steady_clock::now();
    query2(filename);
    toc = std::chrono::steady_clock::now();
    diff = toc - tic;
    std::cout << "Q2," << ncores << "," << diff.count() << std::endl;

    tic = std::chrono::steady_clock::now();
    query3(filename);
    toc = std::chrono::steady_clock::now();
    diff = toc - tic;
    std::cout << "Q3," << ncores << "," << diff.count() << std::endl;

    tic = std::chrono::steady_clock::now();
    query4(filename);
    toc = std::chrono::steady_clock::now();
    diff = toc - tic;
    std::cout << "Q4," << ncores << "," << diff.count() << std::endl;

    tic = std::chrono::steady_clock::now();
    query5(filename);
    toc = std::chrono::steady_clock::now();
    diff = toc - tic;
    std::cout << "Q5," << ncores << "," << diff.count() << std::endl;

    tic = std::chrono::steady_clock::now();
    query6(filename);
    toc = std::chrono::steady_clock::now();
    diff = toc - tic;
    std::cout << "Q6," << ncores << "," << diff.count() << std::endl;

    tic = std::chrono::steady_clock::now();
    query7(filename);
    toc = std::chrono::steady_clock::now();
    diff = toc - tic;
    std::cout << "Q7," << ncores << "," << diff.count() << std::endl;

    tic = std::chrono::steady_clock::now();
    query8(filename);
    toc = std::chrono::steady_clock::now();
    diff = toc - tic;
    std::cout << "Q8," << ncores << "," << diff.count() << std::endl;

    return 0;
}
