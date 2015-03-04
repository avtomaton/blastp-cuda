#!/bin/bash 

#
# Progress Indicator
#
function progress_ind {
  # Sleep at least 1 second otherwise this algoritm will consume
  # a lot system resources.
  interval=1

  while true
  do
    echo -ne "."
    sleep $interval
  done
}

#
# Stop distraction
#
function stop_progress_ind {
  exec 2>/dev/null
  kill $1
  echo -en "\n"
}

ncbi_blast_version="ncbi-blast-2.2.28+-src"
working_dir=`pwd`
nvcc=`which nvcc`

if [ $? -ne 0 ]; then
    echo -e "\nError: nvcc cannot be found. Exiting..."
    exit 1
fi

echo "Do you want to install GPU-BLAST on an existing installation of \"blastp\" [yes/no]"
echo "yes: you will be asked for the installation directory of the \"blastp\" executable"
echo "no: will download and install \"${ncbi_blast_version}\""
read user_input

if [ $user_input == "yes" ]; then
    echo "Please input the installation directory of \"blastp\" of \"${ncbi_blast_version}\""
    read installed_dir
    #check if the executable "blastp" exists in the given directory
    if [ ! -x $installed_dir/blastp ]; then
	echo "Executable \"blastp\" cannot be found"
	echo "Exiting..."
	exit 0
    else
	blastp_version=`$installed_dir/blastp -version | grep blastp | awk '{print $NF}'`
	if [ $blastp_version == "2.2.28+" ]; then 
	    echo "\"blastp\" version $blastp_version is compatible"
	else
	    echo "\"blastp\" version $blastp_version is not supported"
	    echo "Exiting..."
	    exit 0
	fi
    fi
    #check if the source file "blast_engine.c" as a first check that the source code exists
    if [ ! -f "$installed_dir/../../src/algo/blast/core/blast_engine.c" ]; then
	echo "Source code not detected."
	echo "You must have the source code of \"$ncbi_blast_version\" to add GPU-BLAST"
	echo "Exiting..."
	exit 0
    else
	echo "Continuing with the installation of GPU-BLAST..."	
	download=0
    fi
elif [ $user_input == "no" ]; then
    echo "Continuing with the downloading of $ncbi_blast_version..."
    download=1
else
    echo "Please enter \"yes\" or \"no\""
    echo "Exiting.."
    exit 0
fi

if [ $download == 1 ]; then
    
  
#Download the version 2.2.28 of NCBI BLAST
    echo "Downloading NBCI BLAST"
    wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.2.28/${ncbi_blast_version}.tar.gz
    if [ $? -ne 0 ]; then
	echo -e "\nError downloading ${ncbi_blast_version}. Verify that the file exists in the server. Exiting..."
	exit 1
    fi
    echo -e "\nExtracting NCBI BLAST"
    progress_ind &
    pid=$!
    trap "stop_progress_ind $pid; exit" INT TERM EXIT
    
    tar -xzf ${ncbi_blast_version}.tar.gz

    stop_progress_ind $pid

#Move to the c++ directory which contains the configuration script
    cd ${working_dir}/${ncbi_blast_version}/c++/
    
#Configure NCBI BLAST
    progress_ind &
    pid=$!
    trap "stop_progress_ind $pid; exit" INT TERM EXIT
    configuration_options="./configure --without-debug --with-mt --without-sybase --without-ftds --without-fastcgi --without-ncbi-c --without-sssdb --without-sss --without-geo --without-sp --without-orbacus --without-boost"
    echo -e "\nConfiguring ${ncbi_blast_version} with options:"
    echo -e "\n${configuration_options}"
    echo -e "\nIf you want to change these options edit the file \"install\", delete the directory ${ncbi_blast_version}, and rerun install"
    $configuration_options >  ${working_dir}/configure.output
    if [ $? -ne 0 ]; then
	echo -e "\nError in configuring NCBI BLAST installation. Look at \"configure.output\" for more information"
	exit 1
    fi
    
    stop_progress_ind $pid

    #extract the installation directory    
    makefile_line=`cat ../../configure.output | grep Tools`
    makefile=`echo $makefile_line | awk -F" " '{print $NF}'`
    build_dir=`dirname $makefile`
    proj_dir=`dirname $build_dir`
    cxx_dir=`dirname $proj_dir`
else
    build_dir="$installed_dir/../build/"
    proj_dir="$installed_dir/../"
    cxx_dir="$installed_dir/../../"
fi

#Modify NCBI BLAST files
echo -e "\nModifying NCBI BLAST files"
cd ${working_dir}
echo -e "sources in ${cxx_dir}\n"
echo "we are in {$PWD}"
./modify ${cxx_dir} > ${working_dir}/modify.output
echo $?
if [ $? -ne 0 ]; then
    echo -e "\nError: while modirying the source code files. See \"modify.output\" for more details. Exiting..."
    exit 1
fi

#Compile GPU code
progress_ind &
pid=$!
trap "stop_progress_ind $pid; exit" INT TERM EXIT

echo -e "\nCompiling CUDA code"

bin_dir=`dirname ${nvcc}`
cuda_dir=`dirname ${bin_dir}`

os_version=`uname -m`
if [[ "$os_version" = "i686" ]]; then
    cudart_lib=${cuda_dir}/lib
fi

if [[ "$os_version" = "x86_64" ]]; then
    cudart_lib=${cuda_dir}/lib64
fi

if [[ ! -d ${cudart_lib} ]]; then
    echo "CUDA library directory cannot be found. Exiting installation..."
    exit 1
fi

#make -f makefile.shared -B NVCC_PATH=${nvcc} PROJ_DIR=${proj_dir} CXX_DIR=${cxx_dir} CUDA_LIB=${cudart_lib} #> ${working_dir}/gpu_blast.output
echo "project dir: ${proj_dir}"
echo "cxx dir: ${cxx_dir}"
make -f Makefile -B PROJ_DIR=${proj_dir} CXX_DIR=${cxx_dir} > ${working_dir}/gpu_blast.output

#Modify CONF_LIB of Makefile.mk to inlcude the libraries -lgpublast and -lcudart
cp ${build_dir}/Makefile.mk ${build_dir}/Makefile.mk.backup
conf_libs_line=`cat ${build_dir}/Makefile.mk | grep ^CONF_LIB`
conf_libs_line=${conf_libs_line}" -lgpublast -lcudart"
sed "/^CONF_LIBS/ c\ ${conf_libs_line}" ${build_dir}/Makefile.mk > ${build_dir}/Makefile.mk.temp
mv ${build_dir}/Makefile.mk.temp ${build_dir}/Makefile.mk

stop_progress_ind $pid

#Compile NCBI BLAST
echo "Building NCBI BLAST with GPU BLAST"
cd $build_dir
progress_ind &
pid=$!
trap "stop_progress_ind $pid; exit" INT TERM EXIT

make all_r  &>  ${working_dir}/ncbi_blast.output
if [ $? -ne 0 ]; then
    echo -e "\nError in NCBI BLAST installation. Look at \"ncbi_blast.output\" for more information"
    exit 1
fi

stop_progress_ind $pid

exit 0;
