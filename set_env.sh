source ~/envSL7
export SET_ENV_PATH=set_env.sh
#export LD_LIBRARY_PATH=./Delphes${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}
#export ROOT_INCLUDE_PATH=Delphes/external${ROOT_INCLUDE_PATH:+:${ROOT_INCLUDE_PATH}}

#export LD_LIBRARY_PATH=lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}
#export ROOT_INCLUDE_PATH=external_inc${ROOT_INCLUDE_PATH:+:${ROOT_INCLUDE_PATH}}
#export PATH=/usr/local/cuda-11.1/bin:/opt/nvidia/nsight-compute/2020.2.1${PATH:+:${PATH}}

export HEPMC2_INCLUDE_PATH=/homes/jbkim/delphes_madgraph/hepmc2/include
export HEPMC2_LIB_PATH=/homes/jbkim/delphes_madgraph/hepmc2/lib
export HEPMC3_INCLUDE_PATH=/homes/jbkim/delphes_madgraph/hepmc3/include
export HEPMC3_LIB_PATH=/homes/jbkim/delphes_madgraph/hepmc3/lib64
export LD_LIBRARY_PATH=$HEPMC2_LIB_PATH:$HEPMC3_LIB_PATH:lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}
export ROOT_INCLUDE_PATH=external_inc${ROOT_INCLUDE_PATH:+:${ROOT_INCLUDE_PATH}}
export PATH=/usr/local/cuda-11.1/bin:/opt/nvidia/nsight-compute/2020.2.1${PATH:+:${PATH}}

export VIRTUAL_ENV_DISABLE_PROMPT=1
source py-environ/bin/activate
