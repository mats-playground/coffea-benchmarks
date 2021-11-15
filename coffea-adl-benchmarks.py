#!/usr/bin/env python
import os.path
import subprocess
import time
from itertools import product

import awkward as ak
import hist
import numpy as np
import pandas as pd
import psutil
import tqdm

from coffea import processor

# The opendata files are non-standard NanoAOD, so some optional data columns are missing
processor.NanoAODSchema.warn_missing_crossrefs = False

proc = psutil.Process()


def run(query, chunksize, workers, file):
    if not file.startswith("/dev/shm"):
        # https://stackoverflow.com/questions/9551838/how-to-purge-disk-i-o-caches-on-linux
        try:
            subprocess.run("sync", check=True)
            subprocess.run(
                ["sudo", "bash", "-c", "echo 3 > /proc/sys/vm/drop_caches"], check=True
            )
        except PermissionError:
            pass

    tic = time.monotonic()
    cputic = proc.cpu_times()
    output, metrics = processor.run_uproot_job(
        fileset={"SingleMu": [file]},
        treename="Events",
        processor_instance=query(),
        executor=processor.futures_executor
        if workers > 1
        else processor.iterative_executor,
        executor_args={
            "schema": processor.NanoAODSchema,
            "savemetrics": True,
            "status": False,
        },
        chunksize=chunksize,
        maxchunks=None,
    )
    toc = time.monotonic()
    cputoc = proc.cpu_times()
    metrics["query"] = query.__name__
    metrics["tgt_chunksize"] = chunksize
    metrics["chunksize"] = metrics["entries"] / metrics["chunks"]
    metrics["workers"] = workers
    metrics["walltime"] = toc - tic
    metrics["path"] = os.path.dirname(file)
    metrics.update(
        {
            n: f - i
            for n, f, i in zip(
                "user system children_user children_system iowait".split(),
                cputoc,
                cputic,
            )
        }
    )

    return output, metrics


class Q1Processor(processor.ProcessorABC):
    """Plot the <i>E</i><sub>T</sub><sup>miss</sup> of all events."""

    def process(self, events):
        return (
            hist.Hist.new.Reg(100, 0, 200, name="met", label="$E_{T}^{miss}$ [GeV]")
            .Double()
            .fill(events.MET.pt)
        )

    def postprocess(self, accumulator):
        return accumulator


class Q2Processor(processor.ProcessorABC):
    """Plot the <i>p</i><sub>T</sub> of all jets."""

    def process(self, events):
        return (
            hist.Hist.new.Reg(100, 0, 200, name="ptj", label="Jet $p_{T}$ [GeV]")
            .Double()
            .fill(ak.flatten(events.Jet.pt))
        )

    def postprocess(self, accumulator):
        return accumulator


class Q2Kin2DProcessor(processor.ProcessorABC):
    """Plot the <i>p</i><sub>T</sub> of all jets."""

    def process(self, events):
        return (
            hist.Hist.new.Reg(100, 0, 200, name="ptj", label="Jet $p_{T}$ [GeV]")
            .Reg(100, -5, 5, name="etaj", label=r"Jet $\eta$")
            .Double()
            .fill(ak.flatten(events.Jet.pt), ak.flatten(events.Jet.eta))
        )

    def postprocess(self, accumulator):
        return accumulator


class Q2Kin3DProcessor(processor.ProcessorABC):
    """Plot the <i>p</i><sub>T</sub> of all jets."""

    def process(self, events):
        return (
            hist.Hist.new.Reg(100, 0, 200, name="ptj", label="Jet $p_{T}$ [GeV]")
            .Reg(100, -5, 5, name="etaj", label=r"Jet $\eta$")
            .Reg(100, -np.pi, np.pi, name="phij", label=r"Jet $\phi$")
            .Double()
            .fill(
                ak.flatten(events.Jet.pt),
                ak.flatten(events.Jet.eta),
                ak.flatten(events.Jet.phi),
            )
        )

    def postprocess(self, accumulator):
        return accumulator


class Q3Processor(processor.ProcessorABC):
    """Plot the <i>p</i><sub>T</sub> of jets with |<i>η</i>| < 1."""

    def process(self, events):
        return (
            hist.Hist.new.Reg(100, 0, 200, name="ptj", label="Jet $p_{T}$ [GeV]")
            .Double()
            .fill(ak.flatten(events.Jet[abs(events.Jet.eta) < 1].pt))
        )

    def postprocess(self, accumulator):
        return accumulator


class Q4Processor(processor.ProcessorABC):
    """Plot the <i>E</i><sub>T</sub><sup>miss</sup> of events that have at least
    two jets with <i>p</i><sub>T</sub> > 40 GeV.
    """

    def process(self, events):
        has2jets = ak.sum(events.Jet.pt > 40, axis=1) >= 2
        return (
            hist.Hist.new.Reg(100, 0, 200, name="met", label="$E_{T}^{miss}$ [GeV]")
            .Double()
            .fill(events[has2jets].MET.pt)
        )

    def postprocess(self, accumulator):
        return accumulator


class Q5Processor(processor.ProcessorABC):
    """Plot the <i>E</i><sub>T</sub><sup>miss</sup> of events that have an
    opposite-charge muon pair with an invariant mass between 60 and 120 GeV.
    """

    def process(self, events):
        mupair = ak.combinations(events.Muon, 2)
        with np.errstate(invalid="ignore"):
            pairmass = (mupair.slot0 + mupair.slot1).mass
        goodevent = ak.any(
            (pairmass > 60)
            & (pairmass < 120)
            & (mupair.slot0.charge == -mupair.slot1.charge),
            axis=1,
        )
        return (
            hist.Hist.new.Reg(100, 0, 200, name="met", label="$E_{T}^{miss}$ [GeV]")
            .Double()
            .fill(events[goodevent].MET.pt)
        )

    def postprocess(self, accumulator):
        return accumulator


class Q6Processor(processor.ProcessorABC):
    """For events with at least three jets, plot the <i>p</i><sub>T</sub> of the trijet
    four-momentum that has the invariant mass closest to 172.5 GeV in each event and
    plot the maximum <i>b</i>-tagging discriminant value among the jets in this trijet.
    """

    def process(self, events):
        jets = ak.zip(
            {k: getattr(events.Jet, k) for k in ["x", "y", "z", "t", "btag"]},
            with_name="LorentzVector",
            behavior=events.Jet.behavior,
        )
        trijet = ak.combinations(jets, 3, fields=["j1", "j2", "j3"])
        trijet["p4"] = trijet.j1 + trijet.j2 + trijet.j3
        trijet = ak.flatten(
            trijet[ak.singletons(ak.argmin(abs(trijet.p4.mass - 172.5), axis=1))]
        )
        maxBtag = np.maximum(
            trijet.j1.btag,
            np.maximum(
                trijet.j2.btag,
                trijet.j3.btag,
            ),
        )
        return {
            "trijetpt": hist.Hist.new.Reg(
                100, 0, 200, name="pt3j", label="Trijet $p_{T}$ [GeV]"
            )
            .Double()
            .fill(trijet.p4.pt),
            "maxbtag": hist.Hist.new.Reg(
                100, 0, 1, name="btag", label="Max jet b-tag score"
            )
            .Double()
            .fill(maxBtag),
        }

    def postprocess(self, accumulator):
        return accumulator


class Q7Processor(processor.ProcessorABC):
    """Plot the scalar sum in each event of the <i>p</i><sub>T</sub> of jets with
    <i>p</i><sub>T</sub> > 30 GeV that are not within 0.4 in Δ<i>R</i> of any light
    lepton with <i>p</i><sub>T</sub> > 10 GeV.
    """

    def process(self, events):
        cleanjets = events.Jet[
            ak.all(
                events.Jet.metric_table(events.Muon[events.Muon.pt > 10]) >= 0.4, axis=2
            )
            & ak.all(
                events.Jet.metric_table(events.Electron[events.Electron.pt > 10])
                >= 0.4,
                axis=2,
            )
            & (events.Jet.pt > 30)
        ]
        return (
            hist.Hist.new.Reg(
                100, 0, 200, name="sumjetpt", label=r"Jet $\sum p_{T}$ [GeV]"
            )
            .Double()
            .fill(ak.sum(cleanjets.pt, axis=1))
        )

    def postprocess(self, accumulator):
        return accumulator


class Q8Processor(processor.ProcessorABC):
    """For events with at least three light leptons and a same-flavor
    opposite-charge light lepton pair, find such a pair that has the
    invariant mass closest to 91.2 GeV in each event and plot the transverse
    mass of the system consisting of the missing tranverse momentum and
    the highest-<i>p</i><sub>T</sub> light lepton not in this pair.
    """

    def process(self, events):
        events["Electron", "pdgId"] = -11 * events.Electron.charge
        events["Muon", "pdgId"] = -13 * events.Muon.charge
        events["leptons"] = ak.concatenate(
            [events.Electron, events.Muon],
            axis=1,
        )
        events = events[ak.num(events.leptons) >= 3]

        pair = ak.argcombinations(events.leptons, 2, fields=["l1", "l2"])
        pair = pair[(events.leptons[pair.l1].pdgId == -events.leptons[pair.l2].pdgId)]
        with np.errstate(invalid="ignore"):
            pair = pair[
                ak.singletons(
                    ak.argmin(
                        abs(
                            (events.leptons[pair.l1] + events.leptons[pair.l2]).mass
                            - 91.2
                        ),
                        axis=1,
                    )
                )
            ]
        events = events[ak.num(pair) > 0]
        pair = pair[ak.num(pair) > 0][:, 0]

        l3 = ak.local_index(events.leptons)
        l3 = l3[(l3 != pair.l1) & (l3 != pair.l2)]
        l3 = l3[ak.argmax(events.leptons[l3].pt, axis=1, keepdims=True)]
        l3 = events.leptons[l3][:, 0]

        mt = np.sqrt(2 * l3.pt * events.MET.pt * (1 - np.cos(events.MET.delta_phi(l3))))
        return (
            hist.Hist.new.Reg(
                100, 0, 200, name="mt", label=r"$\ell$-MET transverse mass [GeV]"
            )
            .Double()
            .fill(mt)
        )

    def postprocess(self, accumulator):
        return accumulator


if __name__ == "__main__":
    queries = [
        Q1Processor,
        Q2Processor,
        Q3Processor,
        Q4Processor,
        Q5Processor,
        Q6Processor,
        Q7Processor,
        Q8Processor,
    ]
    chunksizes = [2 ** 13, 2 ** 15, 2 ** 17, 2 ** 19, 2 ** 21]
    ncores = [1, 3, 12, 24, 48]
    files = [
        "/dev/shm/Run2012B_SingleMu.root",
        # "/ssd/Run2012B_SingleMu.root",
        # "/magnetic/Run2012B_SingleMu.root",
    ]
    benchpoints = list(product(queries, chunksizes, [24], files))
    benchpoints += list(product(queries, [2**19], ncores, files))
    queries = [Q2Processor, Q2Kin2DProcessor, Q2Kin3DProcessor]
    ncores = [12, 18, 24]
    chunksizes = [2 ** 17, 2 ** 18, 2 ** 19]
    benchpoints += list(product(queries, chunksizes, ncores, files))
    benchpoints = list(set(benchpoints))
    results = []
    for query, chunksize, workers, file in tqdm.tqdm(benchpoints):
        _, metrics = run(query, chunksize, workers, file)
        del metrics["columns"]
        results.append(metrics)

    df = pd.DataFrame(results)
    df["us*core/evt"] = df["walltime"] * 1e6 * df["workers"] / df["entries"]
    df["b/evt"] = df["bytesread"] / df["entries"]
    df["MB/s/core"] = df["bytesread"] * 1e-6 / df["workers"] / df["walltime"]
    print(df)
    df.to_pickle("results.pkl")
