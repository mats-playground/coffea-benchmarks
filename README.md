# Coffea Benchmarks
This repository offers solutions to the benchmarks listed [here](https://github.com/iris-hep/adl-benchmarks-index) completed in [Coffea](https://github.com/CoffeaTeam/coffea).

The file these solutions are intended to run on is public and can be accessed by: `root://eospublic.cern.ch//eos/root-eos/benchmark/Run2012B_SingleMu.root`

As of 1/20/2022, these examples have been contributed by [Nick Smith](https://github.com/nsmith-). They are in the coffea-adl-benchmarks.ipynb notebook and respond to the following queries:

1. Plot the missing ET of all events.
2. Plot pT of all jets in all events.
3. Plot pT of jets with |η| < 1.
4. Plot the missing ET of events that have at least two jets with pT > 40 GeV.
5. Plot the missing ET of events that have an opposite-sign muon pair with an invariant mass between 60 and 120 GeV.
6. Plot pT of the trijet system with the mass closest to 172.5 GeV in each event and plot the maximum b-tagging discriminant value among the jets in the triplet.
7. Plot the sum of pT of jets with pT > 30 GeV that are not within 0.4 in ΔR of any lepton with pT > 10 GeV.
8. For events with at least three leptons and a same-flavor opposite-sign lepton pair, find the same-flavor opposite-sign lepton pair with the mass closest to 91.2 GeV and plot the transverse mass of the missing energy and the leading other lepton.

A more pedagogical approach to solving these examples is contained within the [coffea-casa-tutorials repository](https://github.com/CoffeaTeam/coffea-casa-tutorials/tree/master/examples), intended to be ran on [coffea-casa](https://github.com/CoffeaTeam/coffea-casa).
