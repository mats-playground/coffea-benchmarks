Initial conda install
```bash
sudo yum install tmux git htop iotop
curl -sSL https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -o /tmp/miniconda.sh
sudo bash /tmp/miniconda.sh -bfp /usr/local/
rm -rf /tmp/miniconda.sh
sudo /usr/local/bin/conda update conda
conda create -n bench -c conda-forge conda-build voms jupyterlab dask_labextension python==3.8.* ca-policy-lcg xrootd==5.2.0 "uproot>=4.0.8" coffea==0.7.5 lz4 python-xxhash zstandard git
```

Try to configure voms: copy user `.globus`, configure endpoints:
```bash
$ cat ~/.glite/vomses
"cms" "voms2.cern.ch" "15002" "/DC=ch/DC=cern/OU=computers/CN=voms2.cern.ch" "cms"
"cms" "lcg-voms2.cern.ch" "15002" "/DC=ch/DC=cern/OU=computers/CN=lcg-voms2.cern.ch" "cms"
```
and copied some more stuff into `~/.glite` but was not able to voms-proxy-init with a VO

Copied opendata file
```bash
$ xrdcp root://eospublic.cern.ch//eos/root-eos/benchmark/Run2012B_SingleMu.root /magnetic/
[16.11GB/16.11GB][100%][==================================================][25.46MB/s]
```
and also made a copy in `/ssd/`
