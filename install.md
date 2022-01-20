Initial conda install
```bash
sudo yum install tmux git htop iotop
curl -sSL https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -o /tmp/miniconda.sh
sudo bash /tmp/miniconda.sh -bfp /usr/local/
rm -rf /tmp/miniconda.sh
sudo /usr/local/bin/conda update conda
```

Copied opendata file
```bash
$ xrdcp root://eospublic.cern.ch//eos/root-eos/benchmark/Run2012B_SingleMu.root /magnetic/
[16.11GB/16.11GB][100%][==================================================][25.46MB/s]
```
and also made a copy in `/ssd/` and `/dev/shm`


For coffea benchmarks:
```bash
conda create -n coffea-bench -c conda-forge jupyterlab python==3.8.* ca-policy-lcg coffea==0.7.9 psutil
conda activate coffea-bench
python coffea-adl-benchmarks.py
```

For RDF benchmarks:
```bash
conda create -n root-bench -c conda-forge python==3.8.* root==6.24.6
conda activate root-bench
$CXX -O3 $(root-config --cflags --libs) root_bench.cxx -o root_bench
```
(compilation took 10s)
